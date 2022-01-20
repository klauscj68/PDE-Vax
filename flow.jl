""" 
Compute all the required nonlocal terms for the evolution equations. 
Output is returned as a dictionary of Yℓvℓ's with the integrals stored
in the ylvl's integral entry
"""
function nonlocals!(Y::Solℓvℓ;
		    temp::DSymYℓvℓ=DSymYℓvℓ[:temp=>Yℓvℓ(Y.yˢ.tlvl,undef)],
		    prm::DSymVFl=data())
	nnd = Y.yˢ.tlvl.nnd; t₀ = Y.yˢ.tlvl.t₀[1];
	snds = @view Y.yˢ.tlvl.snds[:];

	# Y nonlocals
	∫line!(Y.yˢ);
	∫line!(Y.yᵛ);
	∫line!(Y.yⁱ);

	# ∫(∂s+∂t)λyˢ w/key ∫∂vλyˢ
	@inbounds for i=1:nnd
		temp[:tempyˢ].ys[i]=∂vλ(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyˢ],Y.yˢ,temp[:∫∂vλyˢ]);
	∫line!(temp[:∫∂vλyˢ]);

	#∫λyˢ
	@inbounds for i=1:nnd
		temp[:tempyˢ].ys[i]=λ(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyˢ],Y.yˢ,temp[:∫λyˢ]);
	∫line!(temp[:∫λyˢ]);

	#∫λ²yˢ
	temp[:tempyˢ].ys[:] = (temp[:tempyˢ].ys).^2;
	∏!(temp[:tempyˢ],Y.yˢ,temp[:∫λ²yˢ]);
	∫line!(temp[:∫λ²yˢ]);

	#∫yˢ
	temp[:∫yˢ].ys[:] = Y.yˢ.ys;
	temp[:∫yˢ].∫yds[1] = Y.yˢ.∫yds[1]

	snds = @view Y.yᵛ.tlvl.snds[:];

	#∫(∂s+∂t)αyᵛ w/key ∫∂vαyᵛ
	@inbounds for i=1:nnd
		temp[:tempyᵛ].ys[i]=∂vα(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyᵛ],Y.yᵛ,temp[:∫∂vαyᵛ]);
	∫line!(temp[:∫∂vαyᵛ]);

	#∫(1-α)αyᵛ,∫(1-α)yᵛ,∫(1-α)²yᵛ w/key ∫Imααyᵛ,∫Imαyᵛ,∫Imα²yᵛ
	@inbounds for i=1:nnd
		temp[:tempyᵛ].ys[i]=α(snds[i],t₀;prm=prm);
		temp[:∫Imααyᵛ].ys[i]=temp[:tempyᵛ].ys[i] |> (x->(1-x)*x);
		temp[:∫Imαyᵛ].ys[i]=1-temp[:tempyᵛ].ys[i];
		temp[:∫Imα²yᵛ].ys[i]=temp[:tempyᵛ].ys[i] |> (x->((1-x)^2));
	end
	∏!(temp[:∫Imααyᵛ],Y.yᵛ,temp[:∫Imααyᵛ]); ∫line!(temp[:∫Imααyᵛ]);
	∏!(temp[:∫Imαyᵛ],Y.yᵛ,temp[:∫Imαyᵛ]); ∫line!(temp[:∫Imαyᵛ]);
	∏!(temp[:∫Imα²yᵛ],Y.yᵛ,temp[:∫Imα²yᵛ]); ∫line!(temp[:∫Imα²yᵛ]);

	snds = @view Y.yⁱ.tlvl.snds[:];

	# ∫βyⁱ
	@inbounds for i=1:nnd
		temp[:tempyⁱ].ys[i] = β(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyⁱ],Y.yⁱ,temp[:∫βyⁱ]);
	∫line!(temp[:∫βyⁱ]);

	#∫(∂s+∂t)βyⁱ w/key ∫∂vβyⁱ
	@inbounds for i=1:nnd
		temp[:tempyⁱ].ys[i]=∂vβ(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyⁱ],Y.yⁱ,temp[:∫∂vβyⁱ]);
	∫line!(temp[:∫∂vβyⁱ]);

	#∫βγyⁱ
	@inbounds for i=1:nnd
		temp[:tempyⁱ].ys[i]=β(snds[i],t₀;prm=prm)*γ(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyⁱ],Y.yⁱ,temp[:∫βγyⁱ]);
	∫line!(temp[:∫βγyⁱ]);
end
function nonlocalsinit!(Y::Solℓvℓ;prm::DSymVFl=data());
	val = nonlocalsinit(;prm=prm);
	nonlocals!(Y;temp=val,prm=prm);
	return val
end

"""
Initialize a dictionary that has keys needed for the nonlocal terms
"""
function nonlocalsinit(;prm::DSymVFl=data())
	nnd = ceil(Int64(prm[:nnd][1]));
	temp = DSymYℓvℓ();

	# yˢ integrals
	nlkeys = [:∫∂vλyˢ,:∫λyˢ,:∫λ²yˢ,:∫yˢ,:tempyˢ];
	tlvl = Tℓvℓ(0.0,nnd,prm[:yˢrg]);
	for key in nlkeys
		temp[key]=Yℓvℓ(tlvl,undef);
	end

	# yᵛ integrals
	nlkeys=[:∫∂vαyᵛ,:∫Imααyᵛ,:∫Imαyᵛ,:∫Imα²yᵛ,:tempyᵛ];
	tlvl = Tℓvℓ(0.0,nnd,prm[:yᵛrg]);
	for key in nlkeys
		temp[key]=Yℓvℓ(tlvl,undef);
	end

	# yⁱ integrals
	nlkeys=[:∫βyⁱ,:∫∂vβyⁱ,:∫βγyⁱ,:tempyⁱ];
	tlvl = Tℓvℓ(0.0,nnd,prm[:yⁱrg]);
	for key in nlkeys
		temp[key]=Yℓvℓ(tlvl,undef);
	end

	return temp
end
function nonlocalsinit(t₀::Float64;prm::DSymVFl=data())
	temp = nonlocalsinit(prm=prm);
	for key in keys(temp)
		temp[key].tlvl.t₀[1] = t₀;
	end
end

"""
Compute the interior flow of the solution at query point
"""
function flow!(s::Float64,t::Float64,
	       Y::Solℓvℓ;
	       prm::DSymVFl=data(),
	       nls::DSymYℓvℓ=nonlocals(Y;prm=prm),
	       temp::VecVw=Vector{Float64}(undef,3))
	# decimal value is 1/√2
	temp[1] = -0.7071067811865475*(λ(s,t;prm=prm)+nls[:∫βyⁱ].∫yds[1])*eval(Y.yˢ,s);
	αval = α(s,t;prm=prm);
	temp[2] = -0.7071067811865475*(αval+(1-αval)*nls[:∫βyⁱ].∫yds[1])*eval(Y.yᵛ,s);
	temp[3] = -0.7071067811865475*γ(s,t;prm=prm)*eval(Y.yⁱ,s);

end
function flow(s::Float64,t::Float64,
              Y::Solℓvℓ;
              prm::DSymVFl=data(),
              nls::DSymYℓvℓ=nonlocals(Y;prm=prm))
	temp = Vector{Float64}(undef,3);
	flow!(s,t,Y;prm=prm,nls=nls,temp=temp);
	return temp
end

"""
Compute the flow for implicit boundary conditions at query point
"""
function ∂flow!(t::Float64,Y::Solℓvℓ;
		prm::DSymVFl=data(),
		nls::DSymYℓvℓ=nonlocals(Y;prm=prm),
		temp::VecVw=Vector{Float64}(undef,3))
	temp[1]=0.0;
	temp[2]=nls[:∫∂vλyˢ].∫yds[1]-nls[:∫λ²yˢ].∫yds[1]-nls[:∫λyˢ].∫yds[1]*nls[:∫βyⁱ].∫yds[1];
	temp[3]=nls[:∫βyⁱ].∫yds[1]*(
			    -nls[:∫λyˢ].∫yds[1]-nls[:∫yˢ].∫yds[1]*nls[:∫βyⁱ].∫yds[1]-nls[:∫∂vαyᵛ].∫yds[1]+
			    (1-α(0.0,t;prm=prm))*Y.yᵛ.ys[1] - nls[:∫Imααyᵛ].∫yds[1] - nls[:∫Imα²yᵛ].∫yds[1]*nls[:∫βyⁱ].∫yds[1]
			   ) + (
				nls[:∫yˢ].∫yds[1]+nls[:∫Imαyᵛ].∫yds[1]
			       )*(
				  nls[:∫∂vβyⁱ].∫yds[1]+β(0.0,t;prm=prm)*Y.yⁱ.ys[1]-nls[:∫βγyⁱ].∫yds[1]
				 );
end
function ∂flow(t::Float64,Y::Solℓvℓ;
               prm::DSymVFl=data(),
               nls::DSymYℓvℓ=nonlocals(Y;prm=prm))
	temp = Vector{Float64}(undef,3);
	∂flow!(t,Y;prm=prm,nls=nls,temp=temp);

	return temp
end

"""
Take an euler step of specified size originating from the given tlvl
"""
function euler!(Δt::Float64,Y::Solℓvℓ;
		temp::Solℓvℓ=deepcopy(Y),
		vtemp::VecVw=Vector{Float64}(undef,3),
		∂vtemp::VecVw=Vector{Float64}(undef,3),
		∂temp::VecVw=Vector{Float64}(undef,2),
		prm::DSymVFl=data(),
		nls::DSymYℓvℓ=nonlocals(Y;prm=prm))
	# Advance the time value
	temp.t₀[1]=Y.t₀[1]+Δt;
	temp.yˢ.tlvl.t₀[1]=Y.yˢ.tlvl.t₀[1]+Δt;
	temp.yᵛ.tlvl.t₀[1]=Y.yᵛ.tlvl.t₀[1]+Δt;
	temp.yⁱ.tlvl.t₀[1]=Y.yⁱ.tlvl.t₀[1]+Δt;

	# Compute the new solution values at the nodals
	#  ∂flow
	∂flow!(Y.t₀[1],Y;prm=prm,nls=nls,temp=∂vtemp);
	temp.yˢ.ys[1]=Y.yˢ.ys[1]+Δt*∂vtemp[1];
	temp.yᵛ.ys[1]=Y.yᵛ.ys[1]+Δt*∂vtemp[2];
	temp.yⁱ.ys[1]=Y.yⁱ.ys[1]+Δt*∂vtemp[3];

	#  flow
	#   accomodate sols w/diff meshes away from s=0 axis
	#   first index is yˢ,yᵛ,yⁱ and second index is node.
	#   assumed all sols have same number of nodes
	gen=[1,1];
	@inbounds for i=1:3*(Y.yˢ.tlvl.nnd-1)
		if gen[2]!=Y.yˢ.tlvl.nnd
			gen[2]+=1;
		else
			gen[1]+=1; gen[2]=2;
		end
		pos=gen[2];comp=gen[1];
		ylvl = ( comp==1 ? Y.yˢ : ( comp==2 ? Y.yᵛ : Y.yⁱ ) );
		templvl = ( comp==1 ? temp.yˢ : ( comp==2 ? temp.yᵛ : temp.yⁱ ) );

		# Find the hyperbolic-∂ point this node originated from and dist
		δχ = hyper∂!(Δt,ylvl.tlvl,ylvl.tlvl.snds[pos],temp.t₀[1];temp=∂temp);
		
		# Compute the flow field at the originating point
		flow!(∂temp[1],∂temp[2],Y;prm=prm,nls=nls,temp=vtemp);

		# Compute the new solution value
		if ∂temp[1]>0.0
			# Solution propagated from the t-level
			templvl.ys[pos] = eval(ylvl,∂temp[1])+δχ*vtemp[comp];
		else
			# Solution propagated from the s-axis
			δt = ∂temp[2]-Y.t₀[1];
			templvl.ys[pos] = ylvl.ys[1]+δt*∂vtemp[comp]+δχ*vtemp[comp];
		end

	end
end

"""
Solve the pde system with implicit, nonlocal ∂-terms
"""
function pdesolve(;prm::DSymVFl=data())
	# Adjust prm to satisfy initial conditions
	data!(prm);

	nsmp = (prm[:Trg][2]-prm[:Trg][1])/prm[:dwnsmp][1] |> ceil |> Int64;
	taxis = convert(Vector{Float64},
			LinRange(prm[:Trg][1],prm[:Trg][2],nsmp));
	δprg = taxis[2]-taxis[1];

	# Initialize the solution at starting time
	ysol = Vector{Solℓvℓ}(undef,nsmp);
	nnd = prm[:nnd][1] |> ceil |> Int64;
	nndsmp = prm[:nndsmp][1] |> ceil |> Int64;
	yˢ = Tℓvℓ(prm[:Trg][1],nnd,prm[:yˢrg]) |> (x->Yℓvℓ(x,undef));
	yᵛ = Tℓvℓ(prm[:Trg][1],nnd,prm[:yᵛrg]) |> (x->Yℓvℓ(x,undef));
	yⁱ = Tℓvℓ(prm[:Trg][1],nnd,prm[:yⁱrg]) |> (x->Yℓvℓ(x,undef));

	@inbounds for i=1:nnd
		yˢ.ys[i]=fˢ(yˢ.tlvl.snds[i];prm=prm);
		yᵛ.ys[i]=fᵛ(yᵛ.tlvl.snds[i];prm=prm);
		yⁱ.ys[i]=fⁱ(yⁱ.tlvl.snds[i];prm=prm);
	end

	sol0  = Solℓvℓ(yˢ,yᵛ,yⁱ); # current level
	sol   = deepcopy(sol0); # single adaptive Euler step
	sol1x = deepcopy(sol0); # first double adaptive euler step
	sol2x = deepcopy(sol0); # double adaptive Euler step	

	# Initialize memory allocations for flow computations
	∂temp = Vector{Float64}(undef,2);
	vtemp = Vector{Float64}(undef,3);
	∂vtemp = Vector{Float64}(undef,3);
	nls = nonlocalsinit!(sol0;prm=prm);
	nls1x = deepcopy(nls);
	
	# Initalize memory allocations for error analysis
	aerr = Vector{Float64}(undef,nnd);
	rerr = Vector{Float64}(undef,nnd);

	# Store solution in ysol vector
	ysol[1] = sol0 |> (x->srefine(nndsmp,x));

	# Integrate the system by an adaptive Euler step
	prg = 0.0; Δt = prm[:dwnsmp][1]; pos=2
	while pos <= nsmp
		flagfd = false;
		while !flagfd
			# single euler step
			euler!(Δt,sol0;
			       temp=sol,∂temp=∂temp,vtemp=vtemp,∂vtemp=∂vtemp,
			       prm=prm,nls=nls);

			# double euler step
			euler!(0.5*Δt,sol0;
			       temp=sol1x,∂temp=∂temp,vtemp=vtemp,∂vtemp=∂vtemp,
			       prm=prm,nls=nls);
			
			nonlocals!(sol1x;temp=nls1x,prm=prm);
			euler!(0.5*Δt,sol1x;
			       temp=sol2x,∂temp=∂temp,vtemp=vtemp,∂vtemp=∂vtemp,
			       prm=prm,nls=nls1x);

			# compute errors and adapt step if accuracy not sufficient
			myerrs!(sol.yˢ.ys,sol2x.yˢ.ys;aerr=aerr,rerr=rerr);
			aerr *= prm[:yˢrg][2]-prm[:yˢrg][1];
			if (maximum(aerr)>prm[:atol][1])&&(maximum(rerr)>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			myerrs!(sol.yᵛ.ys,sol2x.yᵛ.ys;aerr=aerr,rerr=rerr);
			aerr *= prm[:yᵛrg][2]-prm[:yᵛrg][1];
			if (maximum(aerr)>prm[:atol][1])&&(maximum(rerr)>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			myerrs!(sol.yⁱ.ys,sol2x.yⁱ.ys;aerr=aerr,rerr=rerr);
			aerr *= prm[:yⁱrg][2]-prm[:yⁱrg][1];
			if (maximum(aerr)>prm[:atol][1])&&(maximum(rerr)>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			flagfd = true;
		end

		# solution was accepted so prepare for next iteration
		#  store solution if at a downsample
		val = sol2x.t₀[1]-prm[:Trg][1];
		while (val >= prg+δprg)&&(pos<=nsmp)
			ysol[pos] = myinterp(taxis[pos],sol0,sol2x) |> (x->srefine(nndsmp,x));
			pos += 1;
			prg += δprg;
			printprg = prg/(prm[:Trg][2]-prm[:Trg][1]);
			println("Simulation progress: $printprg/1"); 
		end
		sol0 = deepcopy(sol2x); nonlocals!(sol0;temp=nls,prm=prm);
		Δt *= 2;
	end

	return ysol
end

"""
Compute the aggregate recovered population from the output of pdesolve
"""
function yʳsolve!(S::Vector{Solℓvℓ};
		  prm::DSymVFl=data(),
	          temp::VecVw=Vector{Float64}(undef,length(S)),
		  tempyᵛ::Yℓvℓ=deepcopy(S[1].yᵛ),
		  tempyⁱ::Yℓvℓ=deepcopy(S[1].yⁱ))
	n = length(S);
	taxis = [S[i].t₀[1] for i=1:n];

	# Compute derivative vector of yʳ in temp
	@inbounds for i=1:n
		# ∫αyᵛ
		tempyᵛ.ys[:] = α.(S[i].yᵛ.tlvl.snds,taxis[i];prm=prm);
		∏!(tempyᵛ,S[i].yᵛ,tempyᵛ);
		∫line!(tempyᵛ);
		temp[i] = tempyᵛ.∫yds[1];

		# ∫γyⁱ
		tempyⁱ.ys[:] = γ.(S[i].yⁱ.tlvl.snds,taxis[i];prm=prm);
		∏!(tempyⁱ,S[i].yⁱ,tempyⁱ);
		∫line!(tempyⁱ);
		temp[i] += tempyⁱ.∫yds[1];
	end

	# integrate in time by trapezoidal
	fwd = @view temp[2:end];
	prs = @view temp[1:end-1];
	temp[2:end] = 0.5*(fwd+prs)*(taxis[2]-taxis[1]);
	temp[1]=0.0;
	cumsum!(temp,temp);
end
function yʳsolve(S::Vector{Solℓvℓ};
		 prm::DSymVFl=data())
	temp = Vector{Float64}(undef,length(S));
	tempyᵛ = deepcopy(S[1].yᵛ);
	tempyⁱ = deepcopy(S[1].yⁱ)

	yʳsolve!(S;prm=prm,temp=temp,tempyᵛ=tempyᵛ,tempyⁱ=tempyⁱ);

	return temp
end
""" 
Plot the masses of compartments as they evolve in time
"""
function plotm(S::Vector{Solℓvℓ};prm::DSymVFl=data())
	n = length(S);
	@inbounds for i=1:n
		∫line!(S[i].yˢ); ∫line!(S[i].yᵛ); ∫line!(S[i].yⁱ)
	end

	taxis = [S[i].t₀[1] for i=1:n];
	yˢ = [S[i].yˢ.∫yds[1] for i=1:n];
	yᵛ = [S[i].yᵛ.∫yds[1] for i=1:n];
	yⁱ = [S[i].yⁱ.∫yds[1] for i=1:n];

	yʳ = yʳsolve(S;prm=prm);

	Σ = yˢ+yᵛ+yⁱ+yʳ;

	p1 = plot(taxis,[yˢ,yᵛ,yⁱ,yʳ,Σ],labels=["∫yˢds" "∫yᵛds" "∫yⁱds" "∫yʳds" "Σ"],
		  linewidth=3);	
	hline!([1.0+prm[:ρ][1]],linewidth=3,linestyle=:dash,labels="theory Σ")
	plot!(xlabel="time elapsed (days)",ylabel="mass");
	plot!(xtickfontsize=10,ytickfontsize=10,xguidefontsize=12,yguidefontsize=12,size=(400,400));
end
