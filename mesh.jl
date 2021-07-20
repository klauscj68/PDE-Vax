# Library for generating a mesh geometry for solving the vaccination PDE
# conservation law system. Mesh is of form [0,umax] x [0,tmax] in (u,t)
# coordinates
using LinearAlgebra

# node
"""
Structure for storing data associated to a node point in (u,t) space
pt:: coordinates
∂pt:: boundary point originating node along characteristic
δl:: length along characteristic between node point and ∂pt
"""
struct node
	pt::Vector{Float64}
	∂pt::Vector{Float64}
	δl::Float64

	function node(pt::Vector{Float64})
		nnz = sum(pt .< 0);
		@assert nnz == 0 "coordinates must be in quadrant 1"
		
		# Find ∂pt
		pos = pt[1] <= pt[2] ? 1 : 2 # coord det ∂pt
		∂pt = Vector{Float64}(undef,2);
		∂pt[pos] = 0.
		∂pt[pos==1 ? 2 : 1] = maximum(pt) - minimum(pt);

		# Compute δl
		δl = norm(pt-∂pt);

		return new(pt,∂pt,δl)
	end

end

# tri
""" 
Structure for storing triangular elements
ndid:: indices of the triangle nodes relative mesh nd listing
Btr:: barycentric coordinate transformation matrix mapping [1;x;y] into
       [b1;b2;b3] relative the order of nd id
urg:: range of u coordinate values across triangle
trg:: range of t coordinates values across triangle
"""
struct tri
	ndid::Vector{Int64}
	Btr::Matrix{Float64}
	urg::Vector{Float64}
	trg::Vector{Float64}

	function tri(ndid::Vector{Int64},nodes::Vector{node})
		# Define the Bary matrix
		B = [1. 1. 1.;
		     nodes[1].pt nodes[2].pt nodes[3].pt];
		Btr = B^(-1);

		# Compute bounding box for triangle
		max = maximum(B[2:3,:],dims=2);
		min = minimum(B[2:3,:],dims=2);

		# Extract ranges
		urg = [min[1];max[1]];
		trg = [min[2];max[2]];

		# Construct the class element
		return new(ndid,Btr,urg,trg)
	end

end


# mesh
"""
Structure for storing mesh data. A rectangular mesh has each element split
into two triangles along its diagonal which, for us, coincides with 
characteristic
"""
struct mesh
	nd::Vector{node}
	elms::Vector{tri}
	urg::Vector{Float64}
	trg::Vector{Float64}
	ntics::Int64
	saxis::Vector{Float64}
	taxis::Vector{Float64}
	χrg::Vector{Float64}
	χaxis::Vector{Float64}

	function mesh(urg::Vector{Float64},trg::Vector{Float64},ntics::Int64)
		# Compute χrg
		χrg = [0.,sqrt(2)*minimum([urg[2],trg[2]])];

		# Generate s-t axes tics
		saxis = convert(Vector,LinRange(urg[1],urg[2],ntics));
		taxis = convert(Vector,LinRange(trg[1],trg[2],ntics));
		χaxis = convert(Vector,LinRange(χrg[1],χrg[2],ntics));

		## Build vector of nodes
		# Generate node coordinates
		ndu = repeat(reshape(saxis,ntics,1),outer=(1,ntics));
		ndt = repeat(reshape(taxis,1,ntics),outer=(ntics,1));

		ndpts = convert(Matrix,transpose([ndu[:] ndt[:]]));

		# Generate nodes by tic index pairings
		idram = 1:ntics;
		ndidu = repeat(reshape(idram,ntics,1),outer=(1,ntics));
		ndidt = repeat(reshape(idram,1,ntics),outer=(ntics,1));

		ndid = convert(Matrix,transpose([ndidu[:] ndidt[:]]));

		nd = Vector{node}(undef,ntics^2);
		for k=1:ntics^2
			nd[k] = node(ndpts[:,k]);
		end

		## Build vector of elements
		# Build rectanges
		#  Every rectangle is in 1-1 correspondence with bottom left corner
		nrects = (ntics-1)^2
		rects = Matrix{Int64}(undef,4,nrects)
		gen = [0,1];
		for k=1:nrects
			# Grab rectangular indices of bottom left node
			# (matching ndid[:,⋅] in [i;j] form)
			if gen[1] != ntics-1
				gen[1] += 1
			else
				gen[1] = 1;
				gen[2] += 1
			end
			i=gen[1]; j=gen[2];
			
			# Compute (i,j) position in ntics^2 ordering
			# (ie compute ⋅)
			#  q = (j-1)*ntics + i
			#  counter-clockwise tour of rectangle from bottom left
			q = [i,i+1,i+1,i] + [j-1,j-1,j,j]*ntics;

			# Store in rects list
			rects[:,k] = q;
		end
		
		# Build tris for storing elements
		tris = Matrix{tri}(undef,2,nrects);
		for k=1:nrects
			# Build tris associated to this rectangle
			#  Node indices for this triangle
			ram = [rects[1,k],rects[2,k],rects[3,k]];
			t1 = tri(ram,[nd[ram[1]],nd[ram[2]],nd[ram[3]]]);
			
			#  Node indices for second triangle
			ram = [rects[1,k],rects[3,k],rects[4,k]];
			t2 = tri(ram,[nd[ram[1]],nd[ram[2]],nd[ram[3]]]);

			tris[:,k] = [t1,t2];
		end
		elms = tris[:];

		return new(nd,elms,urg,trg,ntics,saxis,taxis,χrg,χaxis)
	end
end

# mbrtri
"""
Test whether given point belongs to given triangle
t:: triangle checking for membership
pt:: point checking to see if in triangle
flagbx:: optional Boolean saying if instead of exact membership
         test to see if can rule out membership based on simple bounding
	 box check
"""
function mbrtri(t::tri,pt::Vector{Float64};flagbx::Bool=false)
	if !flagbx
		# Compute barycentric coordinates
		b = t.Btr*[1.;pt];
		nneg = Float64(sum(b.<0));

		return nneg == 0.
	else
		# Test if belongs to bounding box
		return (t.urg[1] <= pt[1])&(t.urg[1] >= pt[1])&
		        (t.urg[2] <= pt[2])&(t.urg[2] >= pt[2])
	end
end

"""
Find the first triangle containing each query point 
(-1 means not found)
m:: mesh over which testing for membership
q:: 2 x nq array that stores query points
"""
function mbrtri(m::mesh,q::Matrix{Float64})
	nq = size(q)[2];
	nelms = length(m.elms);
	
	# Find triangles containing point
	pos = repeat([-1],outer=(nq,));
	gen = [1,0];
	for k=1:nq*nelms
		# 1<=i<=nq and 1<=j<=nelms
		if gen[2] != nelms
			gen[2] += 1;
		else
			gen[1] += 1;
			gen[2] = 1;
		end
		i = gen[1]; j = gen[2];		
		
		if mbrtri(m.elms[j],q[:,i];flagbx=true)
			pos[i] = mbrtri(m.elms[j],q[:,i]) ? j : -1;
			if pos[i] != -1
				continue
			end
		end
	end
	return pos
end
""" 
Same as other mbrtri but returns an Int64 rather than array since working
with single point
"""
function mbrtri(m::mesh,q::Vector{Float64})
	qram = reshape(q,outer=(2,1));
	pos = mbrtri(m,qram)

	return pos[1]
end

# besttri
"""
Rank triangles closest to the query point by the negativity of their 
barycentric coordinates
"""
function besttri(m::mesh,q::Matrix{Float64})
	nq = size(q)[2];
	nelms = length(m.elms);

	# Find triangles best containing point
	pos = Vector{Int64}(undef,nq);
	gen = [1,0];
	bestval = -Inf;
	bestid = -1;
	for k=1:nq*nelms
		# 1<=i<=nq and 1<=j<=nelms
		if gen[2] != nelms
			gen[2] += 1;
		else
			gen[1] += 1;
			gen[2] = 1;
		end
		i = gen[1]; j = gen[2];
		
		# Reset rankings if a new query point
		if j == 1
			bestval = -Inf;
			bestid = -1;
		end

		# Continue if earlier an exact match was found
		if bestval == Inf
			continue
		end
		
		# Test for membership in triangle
		bq = m.elms[j].Btr*[1.;q[:,i]];
		val = 1. *sum((bq.<0).*bq);
		if val == 0.
			bestid = j;
			pos[i] = bestid;
			bestval = Inf;
			continue
		elseif val > bestval
			bestval = val;
			bestid = j;
		end

		# Store best triangle if finishing with query point
		if j == nelms
			pos[i] = bestid;
		end
	end

	return pos
end
"""
Same as other besttri except because query is a vector, returns a single 
integer rather than an array
"""
function besttri(m::mesh,q::Vector{Float64})
	qram::Matrix{Float64} = reshape(q,(2,1));
	pos = besttri(m,qram);

	return pos[1]
end

