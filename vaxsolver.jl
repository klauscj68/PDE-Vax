## Suite for solving the vaccination PDE system
using DifferentialEquations

#%% Custom Structures
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

	function Domain(srg::Vector{Float64},trg::Vector{Float64},
			nelm::Int64)
		@assert nelm >= 2 "nelm must be at least 2"
		saxis = convert(Vector,LinRange(srg[1],srg[2],nelm));
		taxis = convert(Vector,LinRange(trg[1],trg[2],nelm));

		χaxis = [saxis;taxis[2:end]];

		return new(srg,trg,nelm,saxis,taxis,χaxis)
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
function myfindfirst(tpts::Vector{Float64},teval::Float64)
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
function myinterp(tpts::Vector{Float64},ypts::Vector{Float64},teval::Float64)
	
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
function Fχτ(pt::Union{Vector{Float64},
		       SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}})
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
function Fst(pt::Union{Vector{Float64},
                 SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}})
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
included in quadrature weights for element length. This method sets the default 
nelm value used everywhere throughout the code.
"""
function Gaussb(χτ::Union{Vector{Float64},
		    SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}},
		dom::Domain;
		gaussqd::Dict{Symbol,Vector{Float64}}=quad1d(),
		nelm::Int64=1)
	t = Fχτ(χτ)[2]; L = dom.srg[2];

	npts = length(gaussqd[:b]);
	
	mesh = LinRange(0.,1.,nelm);
	pts = Vector{Matrix{Float64}(2,npts)}(undef,nelm);
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
