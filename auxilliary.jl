## Suite of ancillary routines for solving the vaccination PDE system
# Coordinate system for numerical routines are (s,t)

#%% Custom Structures
# Tℓvℓ
"""
Structure encoding the discretization of [t=t₀] in the (χ,τ)-plane
The inner constructor allows for mutating an optional nds input so
Julia can reduce the memory allocation. t0 is a scalar stored as 
vector so the structures can be mutated
"""
struct Tℓvℓ
	t₀::Vector{Float64}
	nnd::Int64
	srg::Vector{Float64}
	snds::Vector{Float64}
	Δs::Vector{Float64}

	function Tℓvℓ(t₀::Float64,srg::Vector{Float64},nnd::Int64)
		nitv = nnd-1;
		Δs = (srg[end]-srg[1])/(nitv);
		snds = [srg[1]+k*Δs for k=0:nitv];

		return new([t₀],nnd,[srg[1],srg[end]],snds,[Δs])
	end

	function Tℓvℓ(t₀::Float64,nnd::Int64,srg::Vector{Float64})
		return Tℓvℓ(t₀,srg,nnd)
	end

	function Tℓvℓ(x::UndefInitializer,nnd::Int64)
		return new([NaN],nnd,[NaN,NaN],fill(NaN,nnd),[NaN])
	end

	function Tℓvℓ(nnd::Int64,x::UndefInitializer)
		return Tℓvℓ(x,nnd)
	end
end
function Base.show(io::IO,tlvl::Tℓvℓ)
	for fld in fieldnames(Tℓvℓ)
		val = getfield(tlvl,fld);
		flds = string(fld);
		println(io,flds*": $val");
	end
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

	function Yℓvℓ(x::UndefInitializer;nnd::Int64=2)
		tlvl = Tℓvℓ(undef,nnd);
		ys = fill(NaN,nnd);
		∫yds = [NaN];

		return new(tlvl,ys,∫yds)
	end
	function Yℓvℓ(tlvl::Tℓvℓ,ys::Vector{Float64})

		return new(tlvl,ys,[NaN])
	end
	function Yℓvℓ(tlvl::Tℓvℓ,x::UndefInitializer)
		ys = fill(NaN,tlvl.nnd);
		return Yℓvℓ(tlvl,ys)
	end
end
function Base.show(io::IO,ylvl::Yℓvℓ)
	val = ylvl.ys;
	println(io,"ys: $val");
	val = ylvl.∫yds;
	println(io,"∫yds: $val");
	println(io,ylvl.tlvl);
end

# Solℓvℓ
"""
Structure for encoding the solution at the [t=t₀] level
"""
struct Solℓvℓ
	t₀::Vector{Float64}
	yˢ::Yℓvℓ
	yᵛ::Yℓvℓ
	yⁱ::Yℓvℓ
	function Solℓvℓ(yˢ::Yℓvℓ,yᵛ::Yℓvℓ,yⁱ::Yℓvℓ)
		return new(yˢ.tlvl.t₀,yˢ,yᵛ,yⁱ)
	end
	function Solℓvℓ(x::UndefInitializer;nnd::Int64=2)
		yˢ=Yℓvℓ(undef;nnd=nnd);
		yᵛ=Yℓvℓ(undef;nnd=nnd);
		yⁱ=Yℓvℓ(undef;nnd=nnd);

		return Solℓvℓ(yˢ,yᵛ,yⁱ)
	end
end
function Base.show(io::IO,sol::Solℓvℓ)
	val = sol.t₀;
	println(io,"t₀: $val");println("");
	println(io,"yˢ:");println(sol.yˢ);
	println(io,"yᵛ:");println(sol.yᵛ);
	println(io,"yⁱ:");println(sol.yⁱ);
end

#%% Aliases
VecVw = Union{
	      Vector{Float64},
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Int64, Base.Slice{Base.OneTo{Int64}}}, true}, # Slice of Yℓvℓ.Tℓvℓ.nds[i,:]
	      SubArray{Float64, 1, Matrix{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, # Slice of Yℓvℓ.Tℓvℓ.nds[i,:]
	      SubArray{Float64, 1, Vector{Float64}, Tuple{Base.Slice{Base.OneTo{Int64}}}, true}, # Slice of Yℓvℓ.Tℓvℓ.snds[:]
	     };
DSymFl = Dict{Symbol,Float64};
DSymVFl = Dict{Symbol,Vector{Float64}};
DSymYℓvℓ = Dict{Symbol,Yℓvℓ};

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
function myinterp(t::Float64,y1::Yℓvℓ,y2::Yℓvℓ)
	@assert (t>=y1.tlvl.t₀[1])&&(t<=y2.tlvl.t₀[1]) "time levels not sequential in myinterp"
	s = (t-y1.tlvl.t₀[1])/(y2.tlvl.t₀[1]-y1.tlvl.t₀[1]);
	ynew = deepcopy(y1);
	ynew.ys[:] = (1-s)*y1.ys + s*y2.ys;
	ynew.tlvl.t₀[1]=t;

	return ynew
end
function myinterp(t,y1::Solℓvℓ,y2::Solℓvℓ)
	@assert (t>=y1.t₀[1])&&(t<=y2.t₀[1]) "time levels not sequential in myinterp"

	s = (t-y1.t₀[1])/(y2.t₀[1]-y1.t₀[1]);
	ynew = deepcopy(y1);

	ynew.t₀[1] = t;
	ynew.yˢ.tlvl.t₀[1] = t;
	ynew.yᵛ.tlvl.t₀[1] = t;
	ynew.yⁱ.tlvl.t₀[1] = t;

	ynew.yˢ.ys[:] = (1-s)*y1.yˢ.ys + s*y2.yˢ.ys;
	ynew.yᵛ.ys[:] = (1-s)*y1.yᵛ.ys + s*y2.yᵛ.ys;
	ynew.yⁱ.ys[:] = (1-s)*y1.yⁱ.ys + s*y2.yⁱ.ys;

	return ynew
end

# hyper∂
"""
Given the step size along t-axis for which advancing the Euler sol
the tlvl underneath and the query point, find the point on this
defined hyperbolic boundary from which point originated and the 
distance between them
"""
function hyper∂(Δt::Float64,tlvl::Tℓvℓ,qpt::VecVw)
	δt = qpt[1]>Δt ? Δt : qpt[1];
	χpt = qpt .- δt;

	δs = 1.4142135623730951*δt;

	return χpt,δs
end
function hyper∂(tlvl::Tℓvℓ,Δt::Float64,qpt::VecVw)
	return hyper∂(Δt,tlvl,qpt)
end
function hyper∂!(Δt::Float64,tlvl::Tℓvℓ,
		 qpt::VecVw;temp::VecVw=[NaN,NaN])
	δt = qpt[1]>Δt ? Δt : qpt[1];
	temp[:] = qpt .- δt;

	δs = 1.4142135623730951*δt;

	return δs
end
function hyper∂!(tlvl::Tℓvℓ,Δt::Float64,
		 qpt::VecVw;temp::VecVw=[NaN,NaN])
	return hyper∂!(Δt,tlvl,qpt;temp=temp)
end
function hyper∂(Δt::Float64,tlvl::Tℓvℓ,s::Float64,t::Float64)
	δt = s>Δt ? Δt : s;
	χpt = [s,t] .- δt;

	δs = 1.4142135623730951*δt;

	return χpt,δs
end
function hyper∂(tlvl::Tℓvℓ,Δt::Float64,s::Float64,t::Float64)
	return hyper∂(Δt,tlvl,s,t)
end
function hyper∂!(Δt::Float64,tlvl::Tℓvℓ,
		 s::Float64,t::Float64;temp::VecVw=[NaN,NaN])
	δt = s>Δt ? Δt : s;
	temp[1] = s - δt;
	temp[2] = t - δt;

	δs =  1.4142135623730951*δt;

	return δs
end
function hyper∂!(tlvl::Tℓvℓ,Δt::Float64,
		 s::Float64,t::Float64;temp::VecVw=[NaN,NaN])
	return hyper∂!(Δt,tlvl,s,t;temp=temp)
end
function hyper∂!(s::Float64,t::Float64,tlvl::Tℓvℓ,Δt::Float64;
		 temp::VecVw=[NaN,NaN]);
	return hyper∂!(Δt,tlvl,s,t;temp=temp)
end
function hyper∂!(s::Float64,t::Float64,Δt::Float64,tlvl::Tℓvℓ;
		 temp::VecVw=[NaN,NaN]);
	return hyper∂!(Δt,tlvl,s,t;temp=temp)
end

# Solution interpolation
"""
Evaluate the interpolated solution at query points. In first
form, routine interpolates the value at an exact t-lvl using 
just the s-coordinate of the point. In the second form, the 
routine interpolates the solution between levels by interpolating
in between a vector of ylvl's by characteristics (actually that form
needs the vector field so save for later)
"""
function eval(ylvl::Yℓvℓ,s::Float64)
	snds = @view ylvl.tlvl.snds[:];
	yval = @view ylvl.ys[:];

	return myinterp(snds,ylvl.ys,s)
end

# Solution integrations
"""
Integrate the solution at a level by trapezoidal rule
"""
function ∫line(ylvl::Yℓvℓ)
	val = 0.0;
	@inbounds for i=1:(ylvl.tlvl.nnd-1)
		val += 0.5*ylvl.tlvl.Δs[1]*(ylvl.ys[i+1]+ylvl.ys[i]);
	end

	return val
end
function ∫line!(ylvl::Yℓvℓ)
	ylvl.∫yds[1] = ∫line(ylvl);
end

"""
Project a solution onto an alternative space discretization
"""
function srefine(ylvl::Yℓvℓ,nnd::Int64)
	tlvl = Tℓvℓ(ylvl.tlvl.t₀[1],ylvl.tlvl.srg,nnd);
	yvals = [eval(ylvl,s) for s in tlvl.snds];

	return Yℓvℓ(tlvl,yvals);
end
function srefine!(ylvl::Yℓvℓ,nnd::Int64;
		  temp::Yℓvℓ=Yℓvℓ(undef;nnd=nnd))

	temp.tlvl.t₀[1] = ylvl.tlvl.t₀[1];
	temp.tlvl.srg[:] = ylvl.tlvl.srg;
	Δs = (srg[2]-srg[1])/(nnd-1); temp.tlvl.snds[:] = [srg[1]+Δs*k for k=0:nnd-1];
	temp.ylvl.ys[:] = [eval(ylvl,s) for s in temp.tlvl.snds];
	temp.tlvl.Δs[1] = Δs;

end
function srefine(nnd::Int64,ylvl::Yℓvℓ)
	return srefine(ylvl,nnd)
end
function srefine!(nnd::Int64,ylvl::Yℓvℓ;
		  temp::Yℓvℓ=Yℓvℓ(undef;nnd=nnd))
	srefine!(ylvl,nnd;temp=temp)
end
function srefine(nnd::Int64,sol::Solℓvℓ)
	yˢ = srefine(nnd,sol.yˢ);
	yᵛ = srefine(nnd,sol.yᵛ);
	yⁱ = srefine(nnd,sol.yⁱ);

	return Solℓvℓ(yˢ,yᵛ,yⁱ)
end
function srefine(sol::Solℓvℓ,nnd::Int64)
	return srefine(nnd,sol)
end

#%% Arithmetic operators on Yℓvℓ's
"""
Can multiply ylvl's together by multiplying the nodal values. Assumes that
the tlvl's are the same.
"""
function Base.:*(y1::Yℓvℓ,y2::Yℓvℓ)
	return Yℓvℓ(y1.tlvl,y1.ys.*y2.ys)
end
function ∏!(y1::Yℓvℓ,y2::Yℓvℓ,y3::Yℓvℓ)
	y3.ys[:] = y1.ys.*y2.ys;
end


#%% Error analysis
"""
vlow and vhigh are resp meant to be the low and high resolution euler step
solutions. It uses vhigh to compute the relative error
"""
function myerrs!(vlow::VecVw,vhigh::VecVw;
		rerr::VecVw=Vector{Float64}(undef,length(vlow)),
		aerr::VecVw=Vector{Float64}(undef,length(vlow)))
	aerr[:] = abs.(vhigh-vlow);
	rerr[:] = aerr./abs.(vhigh);
end
function myerrs(vlw::VecVw,vhigh::VecVw)
	rerr=Vector{Float64}(undef,length(vlow));
	aerr=Vector{Float64}(undef,length(vlow));

	myerrs!(vlow,vhigh;rerr=rerr,aerr=aerr);

	return aerr,rerr
end
