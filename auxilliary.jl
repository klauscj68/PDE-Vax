## Suite of ancillary routines for solving the vaccination PDE system
#%% Aliases
VecVw = Union{
	      Vector{Float64},
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, # Slice of Yℓvℓ.Tℓvℓ.nds[i,:]
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, # Slice of Yℓvℓ.Tℓvℓ.nds[i,:] 
	     }

#%% Coordinate Transformation
# γℓvℓ
"""
Map parameterizing [t=t₀] over (χ,τ)-plane 
"""
function γℓvℓ(t₀::Float64,χ::Float64)
	 
	pt = [χ;
	     1/√(2)*χ - 1/√(2)*abs(χ) + √(2)*t₀];

	return pt
end

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

#%% Custom Structures
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
	snds::Vector{Float64}
	
	function Tℓvℓ(t₀::Float64,χvals::Vector{Float64};
		      nds::Matrix{Float64}=Matrix{Float64}(undef,2,length(χvals)))
		@assert t₀ >= 0 "t₀ must be nonnegative"
		nnd = length(χvals);
		@assert size(nds) == (2,nnd) "nds must have length as χvals"

		# Generate the nodes and ranges
		snds = Vector{Float64}(undef,nnd);
		@inbounds for i=1:nnd 
			pt = γℓvℓ(t₀,χvals[i]);
			nds[:,i] = pt;
			snds[i] = Fχτ(pt)[1];
		end

		χmin,τmin = minimum(nds,dims=2);
		χmax,τmax = maximum(nds,dims=2);

		χrg = [χmin,χmax];
		τrg = [τmin,τmax];

		return new([t₀],nnd,nds,χrg,τrg,snds);
	end	
end

# Tℓvℓ!
"""
Mutate a given Tℓvℓ structure in the fashion the ODE-Euler integration
scheme will use along characteristics:
Maps the s-coord values to the new [t=t0] level in (χ,τ)-plane
"""
function Tℓvℓ!(t₀::Float64,tlvl::Tℓvℓ)
	# Overwrite fields of tlvl to new values
	tlvl.t₀[1] = t₀;
	
	q = [0.,t₀];
	nds = tlvl.nds;
	snds = tlvl.snds;
	@inbounds for i=tlvl.nnd:-1:2
		q[1] = snds[i];
		nds[:,i] = Fst(q);
	end

	nds[:,1] = [-t₀,0.];

	χmin,τmin = minimum(nds,dims=2);
	χmax,τmax = maximum(nds,dims=2);

	tlvl.χrg[:] = [χmin,χmax];
	tlvl.τrg[:] = [τmin,τmax];
end
function Tℓvℓ!(t₀::Float64,tlvl0::Tℓvℓ,tlvl::Tℓvℓ)
	# Overwrite fields of tlvl to new values
	tlvl.t₀[1] = t₀;
	
	q = [0.,t₀];
	snds0 = tlvl0.snds;
	nds = tlvl.nds;
	tlvl.snds[:] = snds0;
	@inbounds for i=tlvl.nnd:-1:2
		q[1] = snds0[i];
		nds[:,i] = Fst(q);
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
	∫yds::Vector{Float64}

	function Yℓvℓ(x::UndefInitializer)
		tlvl = Tℓvℓ(0.,[0.]);
		ys = [NaN];
		∫yds = [NaN];

		return new(tlvl,ys,∫yds)
	end

	function Yℓvℓ(tlvl::Tℓvℓ,ys::Vector{Float64})

		return new(tlvl,ys,[NaN])
	end
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

"""
A 1d linear interpolation scheme to be used on Yℓvℓ structures. Each
consecutive ypts value is a [t=tpts[i]] level and the interpolation is
essentially done by Δt in the (s,t) plane
NOTE: Assumes the ypts all have same snds values, ie are the same spatial
      discretization
"""
function myinterp(tpts::VecVw,ypts::Vector{Yℓvℓ},teval::Float64)
	if teval <= tpts[1]
		val = deepcopy(ypts[1]);
	elseif teval >= tpts[end]
		val = deepcopy(ypts[end]);
	else
		pos = myfindfirst(tpts,teval);
		t1,t2 = tpts[pos-1:pos];
		s = (teval-t1)/(t2-t1);

		# Construct the Tℓvℓ
		nnd = ypts[1].tlvl.nnd;
		nds = Matrix{Float64}(undef,2,nnd);
		for i=1:nnd
			nds[:,i] = Fst([ypts[1].tlvl.snds[i],teval]);
		end

		χmin,τmin = minimum(nds,dims=2);
		χmax,τmax = maximum(nds,dims=2);

		χrg = [χmin,χmax];
		τrg = [τmin,τmax];

		tlvl = Tℓvℓ(teval,nnd,nds,χrg,τrg,ypts[1].tlvl.snds);

		# Construct the ys
		ys = (1-s)*ypts[pos-1].ys + s*ypts[pos].ys;

		val = Yℓvℓ(tlvl,ys);
	end

	return val
end

# Solution interpolation
function eval(ylvl::Yℓvℓ,χ::Float64)
	χs = @view ylvl.tlvl.nds[1,:];

	return myinterp(χs,ylvl.ys,χ)
end

# Solution integrations
function ∫line(ylvl::Yℓvℓ)
	nds = ylvl.tlvl.nds;
	χs = @view ylvl.tlvl.nds[1,:];
	ys = ylvl.ys;
	nnd = size(nds)[2];

	∫val = 0.;

	@inbounds for i=1:nnd-1
		ndm1 = @view nds[:,i];
		ndp1 = @view nds[:,i+1];
		# Split up integral into [χ<0] and [χ>0] so that in ODE
		# integration δt→0 can make error small ind of space 
		# discretization
		if (χs[i]<0)&&(χs[i+1]>0)
			nd0 = γℓvℓ(ylvl.tlvl.t₀[1],0.);
			f0 = eval(ylvl,0.);

			# Integral in [χ<0]
			Δs = nd0 - ndm1;
			ds = √(Δs[1]^2+Δs[2]^2);
			
			∫val += .5/√(3)*(ys[i]+f0)*ds;

			# Intergal in [χ>=0]
			Δs = ndp1 - nd0;
			ds = √(Δs[1]^2+Δs[2]^2);

			∫val += .5*(f0+ys[i+1])*ds;
		else
			# Interval has same sign
			f1 = ys[i]*(χs[i] >= 0 ? 1 : 1/√(3));
			f2 = ys[i+1]*(χs[i] >= 0 ? 1 : 1/√(3));

			Δs = ndp1 - ndm1;
			ds = √(Δs[1]^2+Δs[2]^2);
			∫val += .5*(f1+f2)*ds;
		end
	end
	
	return ∫val
end
function ∫line!(ylvl::Yℓvℓ)
	∫yds = ∫line(ylvl);
	ylvl.∫yds[1] = ∫yds;
end
