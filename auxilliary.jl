## Suite of ancillary routines for solving the vaccination PDE system

#%% Custom Structures
# Aliases for Views and Vector of views
"""
Convenient alias to accomodate vector or view of vector input
"""
Vw = Union{
	   SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, # matrix slice
	   SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true}, # vector slice
	   SubArray{Float64, 1, Matrix{Float64}, Tuple{UnitRange{Int64}, Int64}, true} # vector slice of matrix slice 
	   };
VecVw = Union{Vector{Float64},Vw}; 
VecVecVw = Union{Vector{Vector{Float64}},
		 Vector{Vw},
		 SubArray{Vector{Float64}, 1, Vector{Vector{Float64}}, Tuple{UnitRange{Int64}}, true}}; # vector of vector slice

# γℓvℓ
"""
Map parameterizing [t=t₀] over (χ,τ)-plane 
"""
function γℓvℓ(t₀::Float64,χ::Float64)
	 
	pt = [χ;
	     1/√(2)*χ - 1/√(2)*abs(χ) + √(2)*t₀];

	return pt
end

# Tℓvℓ
"""
Structure encoding the discretization of [t=t₀] in the (χ,τ)-plane
The inner constructor allows for mutating an optional nds input so
Julia can reduce the memory allocation
"""
struct Tℓvℓ
	t₀::Vector{Float64}
	nnd::Int64
	nds::Matrix{Float64}
	χrg::Vector{Float64}
	τrg::Vector{Float64}
	
	function Tℓvℓ(t₀::Float64,χvals::Vector{Float64};
		      nds::Matrix{Float64}=Matrix{Float64}(undef,2,length(χvals)))
		@assert t₀ >= 0 "t₀ must be nonnegative"
		nnd = length(χvals);
		@assert size(nds) == (2,nnd) "nds must have length as χvals"

		# Generate the nodes and ranges
		@inbounds for i=1:nnd 
			nds[:,i] = γℓvℓ(t₀,χvals[i]);
		end

		χmin,τmin = minimum(nds,dims=2);
		χmax,τmax = maximum(nds,dims=2);

		χrg = [χmin,χmax];
		τrg = [τmin,τmax];

		return new([t₀],nnd,nds,χrg,τrg);
	end	
end

# Tℓvℓ!
"""
Mutate a given Tℓvℓ structure in the fashion the ODE-Euler integration
scheme will use along characteristics:
The first n-1 χ-coordinates of the old nodes become the last n-1 
coordinates of the new nodes and the first node is set at (-t₀,0) in
(χ,τ)-plane
"""
function Tℓvℓ!(t₀::Float64,tlvl::Tℓvℓ)
	# Overwrite fields of tlvl to new values
	tlvl.t₀[1] = t₀;
	
	nds = tlvl.nds;
	@inbounds for i=tlvl.nnd:-1:2
		nds[:,i] = γℓvℓ(t₀,nds[1,i-1]);
	end

	nds[:,1] = [-t₀,0.];

	χmin,τmin = minimum(nds,dims=2);
	χmax,τmax = maximum(nds,dims=2);

	tlvl.χrg[:] = [χmin,χmax];
	tlvl.τrg[:] = [τmin,τmax];
end

# Yℓvℓ
"""
Structure encoding the [t=t₀] level set and the nodal values of the 
for a function defined by the interpolant along the line segment
discretization of the level set
"""
struct Yℓvℓ
	tlvl::Tℓvℓ
	ys::Vector{Float64}
end


#%% Interpolation routines
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

	st = [.5*χ+.5*abs(χ) + 1/√(2)*τ,
	      -.5*χ+.5*abs(χ) + 1/√(2)*τ];

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
	      1/√(2)*s+1/√(2)*t - 1/√(2)*abs(s-t)];

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

# Solution interpolation
function eval(ylvl::Yℓvℓ,χ::Float64)
	χs = @view ylvl.tlvl.nds[1,:];

	return myinterp(χs,ylvl.ys,χ)
end

# Solution integrations
function ∫line(ylvl::Yℓvℓ)
	nds = ylvl.tlvl.nds;
	ys = ylvl.ys;
	nnd = size(nds)[2];

	∫val = 0.;
	@inbounds for i=1:nnd-1
		f1 = ys[i]*(ys[i] >= 0 ? 1 : 1/√(3));
		f2 = ys[i+1]*(ys[i+1] >= 0 ? 1 : 1/√(3));

		Δs = nds[:,i+1] - nds[:,i];
		ds = √(Δs[1]^2+Δs[2]^2);
		∫val += .5*(f1+f2)*ds;
	end
	
	return ∫val
end
