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
		∂pt = pt .- pt[pos];

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

	function mesh(urg::Vector{Float64},trg::Vector{Float64},ntics::Int64)
		# Generate u-t axes tics
		uaxis = LinRange(urg[1],urg[2],ntics); 
		taxis = LinRange(trg[1],trg[2],ntics);

		## Build vector of nodes
		# Generate node coordinates
		ndu = repeat(reshape(uaxis,1,ntics),outer=(ntics,1));
		ndt = repeat(reshape(taxis,ntics,1),outer=(1,ntics));

		ndpts = convert(Matrix,transpose([ndu[:] ndt[:]]));

		# Generate nodes by tic index pairings
		idram = 1:ntics;
		ndidu = repeat(reshape(idram,1,ntics),outer=(ntics,1));
		ndidt = repeat(reshape(idram,ntics,1),outer=(1,ntics));

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
		for k=1:nrects
			# Extract rectangular indices of bottom left node
			# (matching ndid[:,⋅] in [i;j] form)
			#  k = (i-1)*(ntics-1) + j for 1<=i,j<=ntics-1
			#  => k-1 = i¯*(ntics-1)+j¯ for 0<=i¯,j¯<ntics-1
			i = 1 + Int64(floor((k-1)/(ntics-1)));
			j = 1 + (k-1-(i-1)*(ntics-1));
			
			# Compute (i,j) position in ntics^2 ordering
			# (ie compute ⋅)
			#  q = (i-1)*ntics + j
			#  counter-clockwise tour of rectangle from bottom left
			q = [(i-1),i,i,i-1]*ntics + [j,j,j+1,j+1];

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

		return new(nd,elms,urg,trg,ntics)
	end
end
