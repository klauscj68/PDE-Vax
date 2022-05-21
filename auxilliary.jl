using Plots,Measures,Random,CSV,DataFrames,Statistics
gr();
## Suite of ancillary routines for solving the vaccination PDE system
# Coordinate system for numerical routines are (s,t)

#%% Custom Structures
# Tℓvℓ
"""
Structure encoding the discretization of [t=t₀] in the (s,t)-plane.
Structure wraps most fields in vectors for memory allocation
purposes.
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
DStrFl = Dict{String,Float64};
DSymVFl = Dict{Symbol,Vector{Float64}};
DStrVFl = Dict{String,Vector{Float64}};
DSymYℓvℓ = Dict{Symbol,Yℓvℓ};
DSymBool = Dict{Symbol,Bool};

#%% Interpolation routines
# myfindfirst
"""
A binary search routine for finding the first time point greater than the query point
among a given sequence of times. Writing because the findfirst routine is proving to
be expensive in Julia. Routine assumes that tpts is ordered least to great. It returns 
the endpoints when teval falls outside. flagmth indicates whether to use the full binary
routine (can be slow) or whether to use division on an unif grid (can be faster)
"""
function myfindfirst(tpts::VecVw,teval::Float64;
	             flagmth::String="div")
	ntpts = length(tpts);
	
	if teval >= tpts[end];
		return length(tpts);
	elseif teval <= tpts[1]
		return 1;
	end

	# Find the smallest interval of type (,] containing point.
	if flagmth == "div"
		Δt = (tpts[2]-tpts[1]);
		b = (teval-tpts[1])÷Δt + 2 |> Int64;
		return b

	elseif flagmth=="binary"
		a=1; b=ntpts;
		flag_fd = false;

		while !flag_fd
			mid = ceil(Int64,.5*a+.5*b);
			if mid == b
				flag_fd = true;
			elseif teval <= tpts[mid]
				b = mid;
			else
				a = mid;
			end
		end

		return b
	end
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
		s = (teval-tpts[pos-1])/(tpts[pos]-tpts[pos-1]);
		val = ypts[pos-1]+s*(ypts[pos]-ypts[pos-1]);
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
function myinterp(t::Float64,y1::Solℓvℓ,y2::Solℓvℓ)
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
solutions. It uses vhigh to compute the relative error. vgrp's are in place
arguments for storing the solution mean values over their base unit time. 
Routine compute the pointwise average absolute and relative errors in the
solution relative the time units that were specified to them.
"""
function myerrs!(vlow::VecVw,vhigh::VecVw,
		 vgrplow::VecVw,vgrphigh::VecVw;
		 rerr::VecVw=Vector{Float64}(undef,length(vgrplow)),
		 aerr::VecVw=Vector{Float64}(undef,length(vgrplow)),
		 ngrp::Int64=length(vgrplow),
		 ntunit::Int64=(length(vlow)÷ngrp) )	

	vwlow = @view vlow[1:ngrp*ntunit]; vwhigh = @view vhigh[1:ngrp*ntunit];
	rshplow = reshape(vwlow,ngrp,ntunit); rshphigh = reshape(vwhigh,ngrp,ntunit);	
	vgrplow[:] = mean(rshplow,dims=1); vgrphigh[:] = mean(rshphigh,dims=1);

	aerr[:] = abs.(vgrphigh-vgrplow);
	rerr[:] = aerr./abs.(vgrphigh);
end
function myerrs(vlow::VecVw,vhigh::VecVw,
		 vgrplow::VecVw,vgrphigh::VecVw)
	rerr=Vector{Float64}(undef,length(vgrplow));
	aerr=Vector{Float64}(undef,length(vgrplow));

	myerrs!(vlow,vhigh,
		vgrplow,vgrphigh;
		rerr=rerr,aerr=aerr);

	return aerr,rerr
end

"""
Check if error vectors at each node fall within one of the two tolerances.
If no storage is given, it writes to rerr vector.
"""
function myerrtst(aerr::VecVw,rerr::VecVw;
		  prm::DSymVFl=data())	
	@inbounds for i=1:length(aerr)
		if (aerr[i]>prm[:atol][1])&&(rerr[i]>prm[:rtol][1])
			return false
		end
	end
	return true
end

#%% Plotting routine
function RecipesBase.plot(V::Vector{Solℓvℓ})
	# yˢ
	taxis = [Y.t₀[1] for Y in V];
	saxis = V[1].yˢ.tlvl.snds/365;

	z = Matrix{Float64}(undef,length(taxis),length(saxis));
	for i=1:size(z)[1]
		z[i,:] = V[i].yˢ.ys;
	end
	p1 = heatmap(saxis,taxis,z,xlabel="age (years)",ylabel="time elapsed (days)",title="yˢ");

	# yᵛ
	taxis = [Y.t₀[1] for Y in V];
	saxis = V[1].yᵛ.tlvl.snds/31;

	z = Matrix{Float64}(undef,length(taxis),length(saxis));
	for i=1:size(z)[1]
		z[i,:] = V[i].yᵛ.ys;
	end
	p2 = heatmap(saxis,taxis,z,xlabel="time vax'd (months)",ylabel="",title="yᵛ");

	# yⁱ
	taxis = [Y.t₀[1] for Y in V];
	saxis = V[1].yⁱ.tlvl.snds;

	z = Matrix{Float64}(undef,length(taxis),length(saxis));
	for i=1:size(z)[1]
		z[i,:] = V[i].yⁱ.ys;
	end
	p3 = heatmap(saxis,taxis,z,xlabel="time inf'd (days)",ylabel="",title="yⁱ");

	lay = @layout [a b c];
	p = plot(p1,p2,p3,layout=lay,margin=5mm,size=(1200,400))
end
""" 
Plot the masses of compartments as they evolve in time
"""
function RecipesBase.plot(S::Vector{Solℓvℓ},yʳ::VecVw;
		          prm::DSymVFl=data())
	n = length(S);
	@inbounds for i=1:n
		∫line!(S[i].yˢ); ∫line!(S[i].yᵛ); ∫line!(S[i].yⁱ)
	end

	taxis = [S[i].t₀[1] for i=1:n];
	yˢ = [S[i].yˢ.∫yds[1] for i=1:n];
	yᵛ = [S[i].yᵛ.∫yds[1] for i=1:n];
	yⁱ = [S[i].yⁱ.∫yds[1] for i=1:n];

	Σ = yˢ+yᵛ+yⁱ+yʳ;

	p1 = plot(taxis,[yˢ,yᵛ,yⁱ,yʳ,Σ],labels=["∫yˢds" "∫yᵛds" "∫yⁱds" "∫yʳds" "Σ"],
		  linewidth=3);	
	hline!([1.0+prm[:ρ][1]],linewidth=3,linestyle=:dash,labels="theory Σ")
	plot!(xlabel="time elapsed (days)",ylabel="mass");
	plot!(xtickfontsize=10,ytickfontsize=10,xguidefontsize=12,yguidefontsize=12,size=(400,400));
end
"""
Plot the errors between solution and its implicit boundary terms
"""
function plotbd(S::Vector{Solℓvℓ};prm=data())
	data!(prm);	
	n = length(S);
	taxis = [S[i].t₀[1] for i=1:n];

	yˢ = fill(0.0,n);
	yᵛ = similar(yˢ); yⁱ = similar(yˢ);

	# allocate for nonlocals and compute theory ∂-values
	nls = nonlocalsinit(S[1].t₀[1];prm=prm);
	@inbounds for i=1:n
		nonlocals!(S[i];temp=nls,prm=prm);

		yᵛ[i] = nls[:∫λyˢ].∫yds[1];
		yⁱ[i] = (nls[:∫yˢ].∫yds[1]+nls[:∫Imαyᵛ].∫yds[1])*nls[:∫βyⁱ].∫yds[1];
	end

	yˢaerr = abs.([ S[i].yˢ.ys[1]-yˢ[i] for i=1:n]);
	yᵛaerr = abs.([ S[i].yᵛ.ys[1]-yᵛ[i] for i=1:n]);
	yⁱaerr = abs.([ S[i].yⁱ.ys[1]-yⁱ[i] for i=1:n]);

	yᵛ = [ (yᵛ[i] >= prm[:atol][1] ? yᵛ[i] : prm[:atol][1]) for i=1:n];
	yⁱ = [ (yⁱ[i] >= prm[:atol][1] ? yⁱ[i] : prm[:atol][1]) for i=1:n];

	yˢrerr = yˢ;
	yᵛrerr = yᵛaerr./yᵛ;
	yⁱrerr = yⁱaerr./yⁱ;

	yˢaerr *= (prm[:yˢrg][2]-prm[:yˢrg][1]);
	yᵛaerr *= (prm[:yᵛrg][2]-prm[:yᵛrg][1]);
	yⁱaerr *= (prm[:yⁱrg][2]-prm[:yⁱrg][1]);

	p1 = plot(taxis,[yˢaerr,yᵛaerr,yⁱaerr],labels=["yˢ" "yᵛ" "yⁱ"],
		  xlabel="taxis pt",ylabel="error",title="Absolute Errors Implicit Bd",linewidth=3);
	plot!(p1,xguidefontsize=12,yguidefontsize=12,xtickfontsize=10,ytickfontsize=10,
	      legendfontsize=10,margin=5mm);

	p2 = plot(taxis,[yˢrerr,yᵛrerr,yⁱrerr],labels=["yˢ" "yᵛ" "yⁱ"],
		  xlabel="taxis pt",title="Relative Errors Implicit Bd",linewidth=3);
	plot!(p2,xguidefontsize=12,yguidefontsize=12,xtickfontsize=10,ytickfontsize=10,
	      legendfontsize=10,margin=5mm);

	lay = @layout [a b];
	plot(p1,p2,size=(800,400))

end

"""
Plot the given equation coefficient with its various approximations
"""
function RecipesBase.plot(s::Symbol;prm::DSymVFl=data())
	data!(prm);

	if s==:α
		nds = LinRange(prm[:yᵛrg][1],prm[:yᵛrg][2],prm[:nnd][1]|>Int64);
		saxis = LinRange(prm[:yᵛrg][1],prm[:yᵛrg][2],10*prm[:nnd][1]|>Int64);

		αnds = [α(s,0.0;prm=prm) for s in nds];
		αaxis = [α(s,0.0;prm=prm) for s in saxis];

		p₁ = plot(saxis,αaxis,labels="true",linewidth=3,title="α",xlabel="s",ylabel="value");
		plot!(p₁,nds,αnds,labels="discrete",linewidth=3,linestyle=:dash);

		∂αnds = [∂vα(s,0.0;prm=prm) for s in nds];
		∂αaxis = [∂vα(s,0.0;prm=prm) for s in saxis];

		p₂ = plot(saxis,∂αaxis,labels="true",linewidth=3,title="∂α",xlabel="s",ylabel="value");
		plot!(p₂,nds,∂αnds,labels="discrete",linewidth=3,linestyle=:dash);

		lay = @layout [a b];
		plot(p₁,p₂,layout=lay)
	elseif s==:β
		nds = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],prm[:nnd][1]|>Int64);
		saxis = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],10*prm[:nnd][1]|>Int64);

		β₀nds = [β₀(s,0.0;prm=prm) for s in nds];
		βnds = [β(s,0.0;prm=prm) for s in nds];
		β₀axis = [β₀(s,0.0;prm=prm) for s in saxis];
		βaxis = [β(s,0.0;prm=prm) for s in saxis];

		p₁ = plot(saxis,[β₀axis,βaxis],labels=["orig cont" "cap cont"],linewidth=3,title="β",xlabel="s",ylabel="value");
		plot!(p₁,nds,[β₀nds,βnds],labels=["orig discrete" "cap discrete"],linewidth=3,linestyle=:dash);

		∂β₀nds = [∂vβ₀(s,0.0;prm=prm) for s in nds];
		∂βnds = [∂vβ(s,0.0;prm=prm) for s in nds];
		∂β₀axis = [∂vβ₀(s,0.0;prm=prm) for s in saxis];
		∂βaxis = [∂vβ(s,0.0;prm=prm) for s in saxis];

		p₂ = plot(saxis,[∂β₀axis,∂βaxis],labels=["orig cont" "cap cont"],linewidth=3,title="∂β",xlabel="s",ylabel="value");
		plot!(p₂,nds,[∂β₀nds,∂βnds],labels=["orig discrete" "cap discrete"],linewidth=3,linestyle=:dash);

		lay = @layout [a b];
		plot(p₁,p₂,layout=lay)
	elseif s==:γ
		nds = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],prm[:nnd][1]|>Int64);
		saxis = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],10*prm[:nnd][1]|>Int64);

		γ₀nds = [γ₀(s,0.0;prm=prm) for s in nds];
		γnds = [γ(s,0.0;prm=prm) for s in nds];
		γ₀axis = [γ₀(s,0.0;prm=prm) for s in saxis];
		γaxis = [γ(s,0.0;prm=prm) for s in saxis];

		p₁ = plot(saxis,[γ₀axis,γaxis],labels=["orig cont" "cap cont"],linewidth=3,title="γ",xlabel="s",ylabel="value");
		plot!(p₁,nds,[γ₀nds,γnds],labels=["orig discrete" "cap discrete"],linewidth=3,linestyle=:dash);

		∂γ₀nds = [∂vγ₀(s,0.0;prm=prm) for s in nds];
		∂γnds = [∂vγ(s,0.0;prm=prm) for s in nds];
		∂γ₀axis = [∂vγ₀(s,0.0;prm=prm) for s in saxis];
		∂γaxis = [∂vγ(s,0.0;prm=prm) for s in saxis];

		p₂ = plot(saxis,[∂γ₀axis,∂γaxis],labels=["orig cont" "cap cont"],linewidth=3,title="∂γ",xlabel="s",ylabel="value");
		plot!(p₂,nds,[∂γ₀nds,∂γnds],labels=["orig discrete" "cap discrete"],linewidth=3,linestyle=:dash);

		lay = @layout [a b];
		plot(p₁,p₂,layout=lay)
	elseif s==:fˢ
		nds = LinRange(prm[:yˢrg][1],prm[:yˢrg][2],prm[:nnd][1]|>Int64);
		saxis = LinRange(prm[:yˢrg][1],prm[:yˢrg][2],10*prm[:nnd][1]|>Int64);

		fˢnds = [fˢ(s;prm=prm) for s in nds];
		fˢaxis = [fˢ(s;prm=prm) for s in saxis];

		p₁ = plot(saxis,fˢaxis,labels="orig",linewidth=3,title="fˢ",xlabel="s",ylabel="value");
		plot!(p₁,nds,fˢnds,labels="discrete",linewidth=3,linestyle=:dash);
		plot(p₁)
	elseif s==:fᵛ
		nds = LinRange(prm[:yᵛrg][1],prm[:yᵛrg][2],prm[:nnd][1]|>Int64);
		saxis = LinRange(prm[:yᵛrg][1],prm[:yᵛrg][2],10*prm[:nnd][1]|>Int64);

		fᵛnds = [fᵛ(s;prm=prm) for s in nds];
		fᵛaxis = [fᵛ(s;prm=prm) for s in saxis];

		p₁ = plot(saxis,fᵛaxis,labels="orig",linewidth=3,title="fᵛ",xlabel="s",ylabel="value");
		plot!(p₁,nds,fᵛnds,labels="discrete",linewidth=3,linestyle=:dash);
	elseif s==:fⁱ
		nds = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],prm[:nnd][1]|>Int64);
		saxis = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],10*prm[:nnd][1]|>Int64);

		fⁱnds = [fⁱ(s;prm=prm) for s in nds];
		fⁱaxis = [fⁱ(s;prm=prm) for s in saxis];

		p₁ = plot(saxis,fⁱaxis,labels="orig",linewidth=3,title="fⁱ",xlabel="s",ylabel="value");
		plot!(p₁,nds,fⁱnds,labels="discrete",linewidth=3,linestyle=:dash);
	elseif s==:Weibull
		nds = LinRange(prm[:yⁱrg][1],prm[:yⁱrg][2],prm[:nnd][1]|>Int64);
		Wβ = Weibull(prm[:βα][1],prm[:βθ][1]); Wγ = Weibull(prm[:γα][1],prm[:γθ][1]);

		βs = [pdf(Wβ,s) for s in nds]; γs = [pdf(Wγ,s) for s in nds];
		plot(nds,βs,linewidth=3,labels="β");
		plot!(nds,γs,linewidth=3,labels="γ");
		plot!(title="Infectiousness and Recovery Distributions",
		      xlabel="days after infection",ylabel="density")
	end
end

# Routines to save a random number generator
#%% mysaverng
"""
Given a MersenneTwister rng, save its internal states to csv file which may be reloaded for
restarting future runs. Note that rng seed is UInt32 but the CSV will save as Int64. Will
need to convert when reloading.
"""
function mysaverng(rng::MersenneTwister;fname::String="RNG")
	# seed is a Vector{UInt32}
	dftemp = DataFrame("seed"=>rng.seed);
	CSV.write(fname*"seed.csv",dftemp);

	# state is a Random.DSFMT.DSFMT_state which stores a Vector{Int32} at :val
	dftemp = DataFrame("state"=>rng.state.val);
	CSV.write(fname*"state.csv",dftemp);

	# vals is a Vector{Float64}
	dftemp = DataFrame("vals"=>rng.vals);
	CSV.write(fname*"vals.csv",dftemp);

	# ints is a Vector{UInt128}
	dftemp = DataFrame("ints"=>rng.ints);
	CSV.write(fname*"ints.csv",dftemp);

	# idxF is a Int64
	dftemp = DataFrame("idxF"=>[rng.idxF]);
	CSV.write(fname*"idxF.csv",dftemp);

	# idxI is a Int64
	dftemp = DataFrame("idxI"=>[rng.idxI]);
	CSV.write(fname*"idxI.csv",dftemp);

	# adv is a Int64
	dftemp = DataFrame("adv"=>[rng.adv]);
	CSV.write(fname*"adv.csv",dftemp);

	# adv_jump is a BigInt
	dftemp = DataFrame("adv_jump"=>[rng.adv_jump]);
	CSV.write(fname*"adv_jump.csv",dftemp);

	# adv_vals is a Int64
	dftemp = DataFrame("adv_vals"=>[rng.adv_vals]);
	CSV.write(fname*"adv_vals.csv",dftemp);

	# adv_ints is a Int64
	dftemp = DataFrame("adv_ints"=>[rng.adv_ints]);
	CSV.write(fname*"adv_ints.csv",dftemp);

	return
end

#%% myloadrng
"""
Given a MersenneTwister saved into csv's like my mysaverng function, restore a rng with these settings.
fname is the prefix that the CSV files begin with. Note that rng seed is UInt32 but the CSV will save 
as Int64. Will need to convert when reloading.
"""
function myloadrng(;fname::String="RNG")
	# seed is a Vector{UInt32}
	dftemp = CSV.read(fname*"seed.csv",DataFrame);
	myseed = dftemp[!,1] |> (x->convert(Vector{UInt32},x));

	# state is a Random.DSFMT.DSFMT_state which stores a Vector{Int32} at :val
	dftemp = CSV.read(fname*"state.csv",DataFrame);
	myvstate = dftemp[!,1] |> (x->convert(Vector{Int32},x));
	mystate = Random.DSFMT.DSFMT_state(myvstate);

	# vals is a Vector{Float64}
	dftemp = CSV.read(fname*"vals.csv",DataFrame);
	myvals = dftemp[!,1] |> (x->convert(Vector{Float64},x));

	# ints is a Vector{UInt128}
	dftemp = CSV.read(fname*"ints.csv",DataFrame);
	myints = dftemp[!,1] |> (x->convert(Vector{UInt128},x));

	# idxF is a Int64
	dftemp = CSV.read(fname*"idxF.csv",DataFrame);
	myidxF = dftemp[1,1] |> (x->convert(Int64,x));

	# idxI is a Int64
	dftemp = CSV.read(fname*"idxI.csv",DataFrame);
	myidxI = dftemp[1,1] |> (x->convert(Int64,x));

	# adv is a Int64
	dftemp = CSV.read(fname*"adv.csv",DataFrame);
	myadv = dftemp[1,1] |> (x->convert(Int64,x));

	# adv_jump is a BigInt
	dftemp = CSV.read(fname*"adv_jump.csv",DataFrame);
	myadv_jump = dftemp[1,1] |> (x->convert(BigInt,x));

	# adv_vals is a Int64
	dftemp = CSV.read(fname*"adv_vals.csv",DataFrame);
	myadv_vals = dftemp[1,1] |> (x->convert(Int64,x));

	# adv_ints is a Int64
	dftemp = CSV.read(fname*"adv_ints.csv",DataFrame);
	myadv_ints = dftemp[1,1] |> (x->convert(Int64,x));

	return MersenneTwister(myseed,mystate,myvals,myints,myidxF,myidxI,myadv,myadv_jump,myadv_vals,myadv_ints)
end

# Routines to store a dictionary into a vector
# wrtprm
"""
Write the prm dictionary to a column vector for storing to csv's
Uses multiple dispatch
call with no args: returns the dimension of column vector needed to store 
                   and list of keys
call with dictionary etc: returns the column vector stored in order of keys(prm)
"""
function wrtprm()
	prm = data(); vkeys = [key for key in keys(prm)];
	
	# Create a vector of aprp size
	nelm = length(vkeys);

	V = Vector{Float64}(undef,nelm);
	@inbounds for i=1:nelm
		V[i] = prm[vkeys[i]][1];
	end
		       
	return prm,vkeys,V
end
function wrtprm!(prm::Dict{Symbol,Vector{Float64}},vkeys::Vector{Symbol},
                   V::VecVw)

	@inbounds for i=1:length(vkeys)
		V[i] = prm[vkeys[i]][1];
	end
	
end
function wrtprm!(prm1::Dict{Symbol,Vector{Float64}},
		 prm2::Dict{Symbol,Vector{Float64}})
	@inbounds for key in keys(prm1)
		prm2[key][:] = prm1[key];
	end
end

# rdprm
"""
Read a column vector formatted like wrtprm into a dictionary for 
restarting runs
"""
function rdprm(V::Vector{Float64},vkeys::Vector{Symbol})
	prm=Dict{Symbol,Vector{Float64}}();
	@inbounds for i=1:length(vkeys)
		prm[vkeys[i]] = [V[i]];
	end
	
	# Handle keys with multiple components (note these were never varied)
	prm0 = data();
	@inbounds for key in [:yˢrg₀,:yᵛrg₀,:yⁱrg₀,:Trg,:yˢrg,:yᵛrg,:yⁱrg]
		prm[key] = deepcopy(prm0[key]);
	end

	return prm,vkeys
end

#%% Miscellaneous
"""
A normal distribution but rescaled to have height 1 at mean
"""
function mynrm(x::Float64,μ::Float64,σ::Float64)
	return exp(-0.5*(x-μ)^2/σ^2);
end
"""
Derivative of the above
"""
function ∂mynrm(x::Float64,μ::Float64,σ::Float64)
	return mynrm(x,μ,σ)*(-(x-μ)/σ^2)
end
"""
A Heaviside function (greek letter Eta) defined here as
H(x) = 1/2/δ*(x+δ)*1_[|x|≤δ] + 1_[x>δ].
"""
function Ηδ(x::Float64;δ::Float64=1.0)
	if x≤-δ
		return 0.0
	elseif abs(x)≤δ
		return (x+δ)/(2*δ)
	else
		return 1.0
	end
end
"""
A mollified Heaviside function (greek letter Eta) to use for smooth
bump functions. The mollifier is a continuous moving average:
fρ = 1/2/ρ*∫^(x+ρ)_(x-ρ)f(s)ds.
"""
function Ηδρ(x::Float64;δ::Float64=1.0,ρ::Float64=1.0)
	if (x-ρ≤-δ)&&(x+ρ≤-δ)
		return 0.0
	elseif (x-ρ≤-δ)&&(abs(x+ρ)≤δ)
		return (x+ρ+δ)^2/(8*ρ*δ)
	elseif (x-ρ≤-δ)&&(x+ρ>δ)
		return (x+ρ)/(2*ρ)
	elseif (abs(x-ρ)≤δ)&&(abs(x+ρ)≤δ)
		return ( (x+ρ)^2 + 2*δ*(x+ρ) - (x-ρ)^2 - 2*δ*(x-ρ) )/(8*ρ*δ)
	elseif (abs(x-ρ)≤δ)&&(x+ρ>δ)
		return ( 3*δ^2 - (x-ρ)^2 - 2*δ*(x-ρ) )/(8*ρ*δ) + (x+ρ-δ)/(2*ρ)
	else
		return 1.0
	end
end
"""
Derivative of the mollified Heaviside function (greek letter Eta) to use in
bump function calculations
"""
function ∂Ηδρ(x::Float64;δ::Float64=1.0,ρ::Float64=1.0)
	return ( Ηδ(x+ρ;δ=δ)-Ηδ(x-ρ;δ=δ) )/(2*ρ)
end
"""
A mollified bump function over an interval built from the Heaviside functions.
Idea is
1_[a,b] = 1_[a,∞] - 1_[b,∞] ≈ Ηδρ(⋅-a) - Ηδρ(⋅-b)
"""
function ζδρ(x::Float64,
	     a::Float64,b::Float64;δ::Float64=1.0,ρ::Float64=1.0)
	return Ηδρ(x-a;δ=δ,ρ=ρ) - Ηδρ(x-b;δ=δ,ρ=ρ)
end
"""
Derivative of the mollified bump function over an interval [a,b]
built from the Heaviside function
"""
function ∂ζδρ(x::Float64,
	      a::Float64,b::Float64;δ::Float64=1.0,ρ::Float64=1.0)
	return ∂Ηδρ(x-a;δ=δ,ρ=ρ) - ∂Ηδρ(x-b;δ=δ,ρ=ρ)
end
