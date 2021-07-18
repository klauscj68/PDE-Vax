# Collection of scripts for discretizing functions and solutions into
# various spline forms. Requires mesh.jl

# finterp
"""
Structure for interpolating and storing α,β,γ,λ on the rectangular mesh. 
These functions are defined in the inner constructor. The functions are
stored as a vector of nodal values ordered like the nodes of mesh.
"""
struct finterp 
	α::Vector{Float64}
	β::Vector{Float64}
	γ::Vector{Float64}
	λ::Vector{Float64}
	fˢ::Vector{Float64}
	fⁱ::Vector{Float64}

	function finterp(mymesh::mesh)
		nnd = length(mymesh.nd);

		α = Vector{Float64}(undef,nnd);
		β = Vector{Float64}(undef,nnd);
		γ = Vector{Float64}(undef,nnd);
		λ = Vector{Float64}(undef,nnd);
		fˢ= Vector{Float64}(undef,nnd);
		fⁱ= Vector{Float64}(undef,nnd);

		for i=1:nnd
			s,t = mymesh.nd[i];
			α[i] = maximum([(365. -t)/365.,0.]);
			β[i] = .1*s;
			γ[i] = .05*t;
			λ[i] = .05*s*t;
			fˢ[i] = 1.;
			fⁱ[i] = 1.;
		end

		return new(α,β,γ,λ,fˢ,fⁱ)
	end
		
end

# yspl
"""
Structure for storing solution splines via their boundary data
sol:: is either :yˢ,:yᵛ,:yⁱ depending on the case
∂s:: stores nodal values of s-axis boundary values (matching a mesh)
∂t:: stores nodal values of t-axis boudary values (matching a mesh)
"""
struct yspl
	sol::Symbol
	∂s::Vector{Float64}
	∂t::Vector{Float64}
end

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

# ∫fds
"""
Compute the integral of mesh-spline interpolant from saxis[1] to saxis[end]
broken up into elements according to the interior points of saxis.
"""
function ∫fds(c::Vector{Float64},mymesh::mesh,
	      saxis::Vector{Float64},tval::Float64;
	      n::Int64=6,gaussqd::Dict{Symbol,Vector{Float64}}=quad1d(n))
	@assert length(saxis) >= 2 "saxis must define an interval"
	# Evaluate spline at quadrature points
	nelms = length(saxis) -1;
	ceval = Matrix{Float64}(undef,n,nelms);
	for i=1:nelms
		pts = (1. .-gaussqd[:p])*saxis[i]+gaussqd[:p]*saxis[i+1];
		ram = [reshape(pts,(1,n));
		       repeat([tval],outer=(1,n))];

		ceval[:,i] = eval(c,mymesh,ram);
	end
	
	# Evaluate integral
	val = sum(gaussqd[:w].*ceval);

	return val
end
