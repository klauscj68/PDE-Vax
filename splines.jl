# Collection of scripts for discretizing functions and solutions into
# various spline forms. Requires mesh.jl

# quad1d
"""
Output the Guassian quadrature weights and point locations for numerical
integration. n = 6 is exact up to polynomials of degree 11.
n:: number of quadrature points used in the interval
rg:: the interval over which numerically integrating
"""
function quad1d(n::Int64=6)
	if n == 2 
		# trapezoidal rule
		w = [.5,.5];
		p = [0.,1.];
	elseif n == 6
		# Gaussian quadrature
		w = [0.171324492,0.360761573,0.467913935,
		     0.467913935,0.360761573,0.171324492];
		b = [0.033765243,0.169395307,0.380690407,
		     0.619309593,0.830604693,0.966234757];
	else
		@warn "Requested number of quad pts is undefined. "*
		        "Defaulting number ..."
		gaussqd = quad1d();
		w = gaussqd[:w];
		p = gaussqd[:p];
	end

	gaussqd = Dict{Symbol,Vector{Float64}}(:w=>w,:p=>p);
	return gaussqd
end

# myfindfirst
"""
A binary search routine for finding the first time point greater than the query point
among a given sequence of times. Writing because the findfirst routine is proving to
be expensive in Julia. Routine uses that tpts is ordered least to great. It's 
assumed teval is strictly contained in the interval partitioned by tpts.
"""
function myfindfirst(tpts::Array{Float64,1},teval::Float64)
	ntpts = length(tpts);
	
	# Find the smallest interval of type (,] containing point.
	idx = [1,ntpts];
	flag_fd = false;

	while !flag_fd
		mid = ceil(Int64,.5*idx[1]+.5*idx[2]);
		if mid == idx[2]
			flag_fd = true;
		elseif teval <= tpts[mid]
			idx[2] = mid;
		else
			idx[1] = mid;
		end
	end

	return idx[2]
end

#%% myinterp
"""
A simple 1d linear interpolation scheme to extend a discrete data set to an interval

vals: ntpts x 2 array of floats. First column is time, second is function value
teval: time point at which to evaluate
"""
function myinterp(tpts::Array{Float64,1},ypts::Array{Float64,1},teval::Float64)
	
	if teval <= tpts[1]
		val = ypts[1];
	elseif teval >= tpts[end]
		val = ypts[end];
	else
		pos = myfindfirst(tpts,teval);
		t1,t2 = tpts[pos-1:pos];
		s = (teval-t1)/(t2-t1);
		v1 = ypts[pos-1];
		v2 = ypts[pos];
		val = v1+s*(v2-v1);
	end

	return val

end

# eval
"""
Evaluate a 1d linear interpolant spline at given query points
c:: gives nodal coefficient values in the order of the nodes of m
mymesh:: mesh over which evaluating interpolant
q:: 2 x nq array that stores query points
"""
function eval(c::Vector{Float64},mymesh::mesh,q::Matrix{Float64};
	      flagexact::Bool=true)
	nq = size(q)[2];

	# Find triangles containing query points
	pos = flagexact ? mbrtri(mymesh,q) : besttri(mymesh,q);

	# Evaluate spline at queries
	val = Vector{Float64}(undef,nq);
	for k=1:nq
		if pos[k] == -1
			val[k] = NaN;
			continue
		end

		tri = mymesh.elms[pos[k]];
		nds = tri.ndid;
		b = tri.Btr*[1.;q[:,k]];
		cval = [c[nds[1]],c[nds[2]],c[nds[3]]];

		val[k] = sum(cval.*b);
	end

	return val
end
"""
Same as other eval except because working with singple point returns 
Float64 not array
"""
function eval(c::Vector{Float64},mymesh::mesh,q::Vector{Float64})
	qram = reshape(q,outer=(2,1));
	val = eval(c,mymesh,qram);

	return val[1]
end

# ∫fdτ
"""
Compute the 1d integral of mesh-spline interpolant along an unnormalized
direction u. 
c:: nodal values wrt to mymesh of integrand
mymesh:: mesh structure encoding domain geometry
u:: direction along which integrating
uaxis:: element discretization of ∫ᵘᶠᵤᵢspl(pₒ+τu)dτ where uaxis[i],uaxis[i+1]
        defines an element. Within an element the integration is done by 
	quadrature
pₒ:: the base point for directional increment
n:: optional arg, number of quadrature points within an element
gaussqd:: option argument specifying the quadrature rule so that weights dont 
          have to be regenerated each call. Should be consistent with n when
	  supplied
"""
function ∫fdτ(c::Vector{Float64},mymesh::mesh,
	      u::Vector{Float64},uaxis::Vector{Float64},pₒ::Vector{Float64};
	      n::Int64=6,gaussqd::Dict{Symbol,Vector{Float64}}=quad1d(n))
	@assert length(uaxis) >= 2 "uaxis must define an interval"
	# Evaluate spline at quadrature points
	nelms = length(uaxis) -1;
	ceval = Matrix{Float64}(undef,n,nelms);
	for i=1:nelms
		# Virtual 1d points along directional axis of u
		pts = (1. .-gaussqd[:p])*uaxis[i]+gaussqd[:p]*uaxis[i+1];
		# Points in physical space
		ram = pₒ.+u*reshape(pts,(1,n));

		ceval[:,i] = eval(c,mymesh,ram);
	end
	
	# Evaluate integral
	val = sum(gaussqd[:w].*ceval);

	return val
end
