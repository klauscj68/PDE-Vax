## Suite for solving the vaccination PDE system
using DifferentialEquations

#%% Custom Structures
# VecVw
"""
Convenient alias to accomodate vector or view of vector input
"""
VecVw = Union{Vector{Float64},
	       SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}};

# Domain
"""
Structure for declaring the geometry for the initial value problem and its
discretization.
"""
struct Domain
	srg::Vector{Float64} # Interval spanned by saxis
	trg::Vector{Float64} # Interval spanned by taxis
	nelm::Int64 # number of elements along an individual ∂-axis
	saxis::Vector{Float64} # Location of saxis nodes along srg
	taxis::Vector{Float64} # Location of taxis nodes along trg
	χaxis::Vector{Float64} # (mirrored taxis) and saxis combined
	nnd::Int64

	function Domain(srg::Vector{Float64},trg::Vector{Float64},
			nelm::Int64)
		@assert nelm >= 2 "nelm must be at least 2"
		saxis = convert(Vector,LinRange(srg[1],srg[2],nelm));
		taxis = convert(Vector,LinRange(trg[1],trg[2],nelm));

		χaxis = [saxis;taxis[2:end]];
		nnd = length(χaxis);

		return new(srg,trg,nelm,saxis,taxis,χaxis,nnd)
	end
end

#%% Ancillary routines
# quad1d
"""
Output the Guassian quadrature weights and point locations for numerical
integration. n = 6 is exact up to polynomials of degree 11. In all cases,
the weights should be multiplied by length of interval over which you are 
integrating.
n:: number of quadrature points used in the interval
"""
function quad1d(n::Int64=6)
	if n == 2 
		# trapezoidal rule
		w = [.5,.5];
		b = [0.,1.];
	elseif n == 6
		# Gaussian quadrature
		w = [0.171324492,0.360761573,0.467913935,
		     0.467913935,0.360761573,0.171324492];
		w *= .5;
		b = [0.033765243,0.169395307,0.380690407,
		     0.619309593,0.830604693,0.966234757];
	else
		@warn "Requested number of quad pts is undefined. "*
		        "Defaulting number ..."
		gaussqd = quad1d();
		w = gaussqd[:w];
		b = gaussqd[:b];
	end

	gaussqd = Dict{Symbol,Vector{Float64}}(:w=>w,:b=>b);
	return gaussqd
end

# myfindfirst
"""
A binary search routine for finding the first time point greater than the query point
among a given sequence of times. Writing because the findfirst routine is proving to
be expensive in Julia. Routine assumes that tpts is ordered least to great. It returns 
the endpoints when teval falls outside.
"""
function myfindfirst(tpts::VecVw,teval::Float64)
	ntpts = length(tpts);
	
	if teval >= tpts[end];
		return length(tpts);
	elseif teval <= tpts[1]
		return 1;
	end

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

# myinterp
"""
A simple 1d linear interpolation scheme to extend a discrete data set to an interval

vals: ntpts x 2 array of floats. First column is time, second is function value
teval: time point at which to evaluate
"""
function myinterp(tpts::VecVw,ypts::VecVw,teval::Float64)
	
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

#%% Coordinate Transformation
# Fχτ
"""
Map from (χ,τ) coordinates to (s,t) coordinates
"""
function Fχτ(pt::VecVw)
	χ = pt[1]; τ = pt[2];

	st = [.5*χ+.5*abs(χ) + 1/sqrt(2)*τ,
	      -.5*χ+.5*abs(χ) + 1/sqrt(2)*τ];

	return st
end
function Fχτ(pts::Matrix{Float64})
	@assert size(pts)[1] == 2 "pts must be two dimensional"
	npts = size(pts)[2];
	
	sts = Matrix{Float64}(undef,2,npts)
	@inbounds for j=1:npts
		sts[:,j] = Fχτ(@view pts[:,j]);
	end

	return sts
end

# Fst
"""
Map from (s,t) coordinates to (χ,τ) coordinates
"""
function Fst(pt::VecVw)
	s = pt[1]; t = pt[2];

	χτ = [s - t,
	      1/sqrt(2)*s+1/sqrt(2)*t - 1/sqrt(2)*abs(s-t)];

	return χτ
end
function Fst(pts::Matrix{Float64})
	@assert size(pts)[1] == 2 "pts must be two dimensional"
	npts = size(pts)[2];
	
	χτs = Matrix{Float64}(undef,2,npts)
	@inbounds for j=1:npts
		χτs[:,j] = Fst(@view pts[:,j]);
	end

	return χτs
end

# Gaussb
"""
Map the gaussian quadrature points in [0,1] to their evaluation point in 
the (χ,τ) plane. These terms come from the pullback of the righthand side
of the system to (χ,τ) coordinates and accordingly depend on those parameters.
Returns a Vector{Matrix{Float64}} where [i][:,j] entry records location of
jᵗʰ quadrature point within iᵗʰ element.

Note: if nelm != 1, then extra Jacobian factor (x 1/nelm) will need to be
included in quadrature weights for element length. 
"""
function Gaussb(χτ::VecVw,
		dom::Domain;
		gaussqd::Dict{Symbol,Vector{Float64}}=quad1d(),
		nelm::Int64)
	t = Fχτ(χτ)[2]; L = dom.srg[2];

	npts = length(gaussqd[:b]);
	
	mesh = LinRange(0.,1.,nelm);
	pts = Vector{Matrix{Float64}}(undef,nelm);
	@inbounds for i=1:nelm
		pts[i] = Matrix{Float64}(undef,2,npts);
	end
	
	gen = [1,0];
	@inbounds for k=1:npts*nelm
		# cycle generator
		if gen[2] != npts
			gen[2] += 1
		else
			gen[1] += 1;
			gen[2] = 1;
		end
		i = gen[1]; j = gen[2];	
		
		# Map barycentric coordinates into element
		ν = mesh[i]*(1-gaussqd[:b][j]) + mesh[i+1]*gaussqd[:b][j];
		
		χeff = -t + L*ν;
		pts[i][:,j] = [χeff,
			       1/sqrt(2)*χeff - 1/sqrt(2)*abs(χeff)+sqrt(2)*t];
	end

	return pts
end
function Gaussb(χτ::Matrix{Float64},dom::Domain;
		gaussqd::Dict{Symbol,Vector{Float64}}=quad1d())
	@assert size(χτ)[1] == 2 "pts must be two dimensional"
	npts = size(χτ)[2];

	bpts = Vector{Vector{Matrix{Float64}}}(undef,npts);
	@inbounds for j=1:npts
		bpts[j] = Gaussb(@view χτ[:,j],dom;gaussqd=gaussqd);
	end

	return bpts
end

#%% Evaluation
"""
Given nodal values c over nodes dom.χaxis in [-T,L], evaluate interpolant
spline at query point χ.  You are thinking of c(τ) already evaluated for a
specific τ and now you are left to evaluate the space component.
"""
function eval(c::VecVw,χ::Float64,dom::Domain)
	@assert length(c) == dom.nnd "length of c does not match nodal dof"
	val = myinterp(dom.χaxis,c,χ);

	return val
end
function eval(c::VecVw,χs::VecVw,dom::Domain)
	npts = length(χs);
	
	vals = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		vals[i] = eval(c,χs[i],dom);
	end

	return vals
end
function eval(c::VecVw,dom::Domain,χ::Float64)
	val = eval(c,χ,dom);
end
function eval(c::VecVw,dom::Domain,χs::VecVw)
	vals = eval(c,χs,dom);
end

#%% Line Integration
# pullb∫fds!
"""
Given a forcing term encoded as a time series of nodal values f in the 
(χ,τ) coordinate plane compute its corresponding pullback line integral:
∫ᴸ₀f(s,t)ds 
at the pt χτ = (χ,τ). The key is to write this integral in terms of f\circ Fχτ
since these are the quantities described by the change of coordinates ODE 
sytem.

f:: Vector of time series nodal values along dom.χaxis
τs:: Series of τ-values for which f values correspond
χτ:: gives χ,τ values at which pulling back integral
gaussqd:: output of quad1d for computing gaussian quadrature
nelm:: Number of subelements for discretizing the [0,1] reference element ∫
       This method sets the default nelm value used everywhere throughout 
       the code.
ramGb:: Optional input that can be used to preallocate an array that is rewritten
        at every iteration of the integration and to not have this memory be
	reallocated each time routine is called. This is only mutated argument.
"""
function pullb∫fds!(f::Vector{VecVw},τs::VecVw,pt::VecVw,
		   dom::Domain;
		   gaussqd = quad1d(),
		   nelm::Int64=1,
		   ramGb::Vector{Matrix{Float64}}=Vector{Matrix{Float64}}(undef,nelm))
	χ = pt[1]; τ = pt[2]; δL = dom.saxis[2]/nelm;	
	
	# Integrate over element divisions of the [0,1] reference element
	∫f = 0. 
	for i=1:nelm
		# extract the gaussian quadrature points  
		
end
