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

	#∫(1-α)αyᵛ,∫(1-α)yᵛ,∫(1-α)²yᵛ w/key ∫Imααyᵛ,∫Imαyᵛ,∫Imα²yᵛ,∫αyᵛ
	@inbounds for i=1:nnd
		temp[:tempyᵛ].ys[i]=α(snds[i],t₀;prm=prm);
		temp[:∫Imααyᵛ].ys[i]=temp[:tempyᵛ].ys[i] |> (x->(1-x)*x);
		temp[:∫Imαyᵛ].ys[i]=1-temp[:tempyᵛ].ys[i];
		temp[:∫Imα²yᵛ].ys[i]=temp[:tempyᵛ].ys[i] |> (x->((1-x)^2));
	end
	∏!(temp[:∫Imααyᵛ],Y.yᵛ,temp[:∫Imααyᵛ]); ∫line!(temp[:∫Imααyᵛ]);
	∏!(temp[:∫Imαyᵛ],Y.yᵛ,temp[:∫Imαyᵛ]); ∫line!(temp[:∫Imαyᵛ]);
	∏!(temp[:∫Imα²yᵛ],Y.yᵛ,temp[:∫Imα²yᵛ]); ∫line!(temp[:∫Imα²yᵛ]);
	∏!(temp[:tempyᵛ],Y.yᵛ,temp[:∫αyᵛ]); ∫line!(temp[:∫αyᵛ]);

	snds = @view Y.yⁱ.tlvl.snds[:];

	#∫(∂s+∂t)βyⁱ w/key ∫∂vβyⁱ
	@inbounds for i=1:nnd
		temp[:tempyⁱ].ys[i]=∂vβ(snds[i],t₀;prm=prm);
	end
	∏!(temp[:tempyⁱ],Y.yⁱ,temp[:∫∂vβyⁱ]);
	∫line!(temp[:∫∂vβyⁱ]);

	# ∫β (intermediate)
	@inbounds for i=1:nnd
		temp[:∫βyⁱ].ys[i] = β(snds[i],t₀;prm=prm);
	end

	# ∫γ (intermediate)
	@inbounds for i=1:nnd
		temp[:∫γyⁱ].ys[i]=γ(snds[i],t₀;prm=prm);
	end

	# ∫βγyⁱ
	∏!(temp[:∫βyⁱ],temp[:∫γyⁱ],temp[:tempyⁱ]);
	∏!(temp[:tempyⁱ],Y.yⁱ,temp[:∫βγyⁱ]);
	∫line!(temp[:∫βγyⁱ]);

	# ∫βyⁱ
	∏!(temp[:∫βyⁱ],Y.yⁱ,temp[:∫βyⁱ]);
	∫line!(temp[:∫βyⁱ]);

	# ∫γyⁱ
	∏!(temp[:∫γyⁱ],Y.yⁱ,temp[:∫γyⁱ]);
	∫line!(temp[:∫γyⁱ]);

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
	nlkeys=[:∫∂vαyᵛ,:∫Imααyᵛ,:∫Imαyᵛ,:∫Imα²yᵛ,:∫αyᵛ,:tempyᵛ];
	tlvl = Tℓvℓ(0.0,nnd,prm[:yᵛrg]);
	for key in nlkeys
		temp[key]=Yℓvℓ(tlvl,undef);
	end

	# yⁱ integrals
	nlkeys=[:∫βyⁱ,:∫∂vβyⁱ,:∫βγyⁱ,:∫γyⁱ,:tempyⁱ];
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
	return temp
end

"""
Compute the interior flow of the solution at query point
"""
function flow!(s::Float64,t::Float64,
	       Y::Solℓvℓ;
	       prm::DSymVFl=data(),
	       nls::DSymYℓvℓ=nonlocals(Y;prm=prm),
	       temp::VecVw=Vector{Float64}(undef,3),
	       comp::Int64=1)
	# decimal value is 1/√2
	if comp==1
		temp[1] = -0.7071067811865475*(λ(s,t;prm=prm)+nls[:∫βyⁱ].∫yds[1])*eval(Y.yˢ,s);
	elseif comp==2
		αval = α(s,t;prm=prm);
		temp[2] = -0.7071067811865475*(αval+(1-αval)*nls[:∫βyⁱ].∫yds[1])*eval(Y.yᵛ,s);
	elseif comp==3
		temp[3] = -0.7071067811865475*γ(s,t;prm=prm)*eval(Y.yⁱ,s);
	end
end
function flow(s::Float64,t::Float64,
              Y::Solℓvℓ;
              prm::DSymVFl=data(),
              nls::DSymYℓvℓ=nonlocals(Y;prm=prm),
	      comp::Int64=1)
	temp = Vector{Float64}(undef,3);
	flow!(s,t,Y;prm=prm,nls=nls,temp=temp,comp=comp);
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
Compute the yʳ flow 
"""
function yʳflow!(t::Float64;
		 prm::DSymVFl=data(),
		 nls::DSymYℓvℓ=nonlocals(Y;prm=prm),
		 temp::VecVw=[0.0])
	temp[:] = nls[:∫αyᵛ].∫yds + nls[:∫γyⁱ].∫yds;
end
function yʳflow(t::Float64;
		prm::DSymVFl=data(),
		nls::DSymYℓvℓ=nonlocals(Y;prm=prm))
	temp = [0.0]
	yʳflow(t;prm=prm,nls=nls,temp=temp);

	return temp
end
"""
Take an euler step of specified size originating from the given tlvl
"""
function euler!(Δt::Float64,Y::Solℓvℓ,yʳ::VecVw;
		temp::Solℓvℓ=deepcopy(Y),
		vtemp::VecVw=Vector{Float64}(undef,3),
		∂vtemp::VecVw=Vector{Float64}(undef,3),
		∂temp::VecVw=Vector{Float64}(undef,2),
		yʳtemp::VecVw=Vector{Float64}(undef,1),
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
	gen=[0,2];
	@inbounds for i=1:3*(Y.yˢ.tlvl.nnd-1)
		if gen[1]!=3
			gen[1]+=1; 
		else
			gen[1]=1; gen[2]+=1;
		end	
		pos=gen[2];comp=gen[1];
		ylvl = ( comp==1 ? Y.yˢ : ( comp==2 ? Y.yᵛ : Y.yⁱ ) );
		templvl = ( comp==1 ? temp.yˢ : ( comp==2 ? temp.yᵛ : temp.yⁱ ) );

		# Find the hyperbolic-∂ point this node originated from and dist
		δχ = hyper∂!(Δt,ylvl.tlvl,ylvl.tlvl.snds[pos],temp.t₀[1];temp=∂temp);
		
		# Compute the flow field at the originating point
		flow!(∂temp[1],∂temp[2],Y;prm=prm,nls=nls,temp=vtemp,comp=comp);

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

	#  yʳ flow
	yʳflow!(Y.t₀[1];prm=prm,nls=nls,temp=yʳtemp);
	yʳtemp[:] = yʳ+Δt*yʳtemp;
end

"""
Solve the pde system with implicit, nonlocal ∂-terms. The flag
says whether progress through the pde should be printed to stdout
"""
function pdesolve(;prm::DSymVFl=data(),
	           flagprg::Bool=true)
	# Adjust prm to satisfy initial conditions
	data!(prm);

	nsmp = (prm[:Trg][2]-prm[:Trg][1])/prm[:dwnsmp][1] |> ceil |> Int64;
	taxis = convert(Vector{Float64},
			LinRange(prm[:Trg][1],prm[:Trg][2],nsmp));
	δprg = taxis[2]-taxis[1];

	# Initialize the solution at starting time
	ysol = Vector{Solℓvℓ}(undef,nsmp);
	yʳsol = Vector{Float64}(undef,nsmp);
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

	yʳ0 = [0.0];
	yʳ = [0.0];
	yʳ1x = [0.0];
	yʳ2x = [0.0];

	# Initialize memory allocations for flow computations
	∂temp = Vector{Float64}(undef,2);
	vtemp = Vector{Float64}(undef,3);
	∂vtemp = Vector{Float64}(undef,3);
	yʳtemp = Vector{Float64}(undef,1);
	nls = nonlocalsinit!(sol0;prm=prm);
	nls1x = deepcopy(nls); nls2x = deepcopy(nls);
	
	# Initalize memory allocations for error analysis
	#  The errors look at error in the avg per age year/per day for susceptible, vax'd and inf	
	#   Note: ngrp means number nodes w/in a tunit
	yˢntunit = (prm[:yˢrg][2]-prm[:yˢrg][1])÷365 |> Int64; yˢngrp = nnd÷yˢntunit; 
	yᵛntunit = (prm[:yᵛrg][2]-prm[:yᵛrg][1])÷1 |> Int64; yᵛngrp = nnd÷yᵛntunit;
	yⁱntunit = (prm[:yⁱrg][2]-prm[:yⁱrg][1])÷1 |> Int64; yⁱngrp = nnd÷yⁱntunit;
	
	#aerr = Vector{Float64}(undef,nnd); rerr = Vector{Float64}(undef,nnd);
	yˢaerr = Vector{Float64}(undef,yˢntunit); yˢrerr = Vector{Float64}(undef,yˢntunit);
	yˢlow = Vector{Float64}(undef,yˢntunit); yˢhigh = Vector{Float64}(undef,yˢntunit);

	yᵛaerr = Vector{Float64}(undef,yᵛntunit); yᵛrerr = Vector{Float64}(undef,yᵛntunit);
	yᵛlow = Vector{Float64}(undef,yᵛntunit); yᵛhigh = Vector{Float64}(undef,yᵛntunit);

	yⁱaerr = Vector{Float64}(undef,yⁱntunit); yⁱrerr = Vector{Float64}(undef,yⁱntunit);
	yⁱlow = Vector{Float64}(undef,yⁱntunit); yⁱhigh = Vector{Float64}(undef,yⁱntunit);

	# Store solution in ysol vector
	ysol[1] = sol0 |> (x->srefine(nndsmp,x));
	yʳ[1] = 0.0;

	# Integrate the system by an adaptive Euler step
	prg = 0.0; Δt = prm[:dwnsmp][1]; pos=2; nΔtfail = 0;
	while pos <= nsmp
		flagfd = false;
		while !flagfd	

			# single euler step
			euler!(Δt,sol0,yʳ0;
			       temp=sol,∂temp=∂temp,vtemp=vtemp,∂vtemp=∂vtemp,yʳtemp=yʳ,
			       prm=prm,nls=nls);

			# double euler step
			euler!(0.5*Δt,sol0,yʳ0;
			       temp=sol1x,∂temp=∂temp,vtemp=vtemp,∂vtemp=∂vtemp,yʳtemp=yʳ1x,
			       prm=prm,nls=nls);
			
			nonlocals!(sol1x;temp=nls1x,prm=prm);
			euler!(0.5*Δt,sol1x,yʳ1x;
			       temp=sol2x,∂temp=∂temp,vtemp=vtemp,∂vtemp=∂vtemp,yʳtemp=yʳ2x,
			       prm=prm,nls=nls1x);

			nonlocals!(sol2x;temp=nls2x,prm=prm);

			if Δt<=prm[:Δtmin][1]
				nΔtfail += 1;
				break
			end

			#∫line!(sol.yˢ); ∫line!(sol.yᵛ); ∫line!(sol.yⁱ);
			#∫line!(sol2x.yˢ); ∫line!(sol2x.yᵛ); ∫line!(sol2x.yⁱ);
			# compute errors and adapt step if accuracy not sufficient
			#  nonlocal errors
			aerr = abs(nls2x[:∫∂vλyˢ].∫yds[1]-nls[:∫∂vλyˢ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫∂vλyˢ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫λyˢ].∫yds[1]-nls[:∫λyˢ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫λyˢ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫λ²yˢ].∫yds[1]-nls[:∫λ²yˢ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫λ²yˢ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫yˢ].∫yds[1]-nls[:∫yˢ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫yˢ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫∂vαyᵛ].∫yds[1]-nls[:∫∂vαyᵛ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫∂vαyᵛ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫Imααyᵛ].∫yds[1]-nls[:∫Imααyᵛ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫Imααyᵛ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫Imαyᵛ].∫yds[1]-nls[:∫Imαyᵛ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫Imαyᵛ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫Imα²yᵛ].∫yds[1]-nls[:∫Imα²yᵛ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫Imα²yᵛ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫αyᵛ].∫yds[1]-nls[:∫αyᵛ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫αyᵛ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫βyⁱ].∫yds[1]-nls[:∫βyⁱ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫βyⁱ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫∂vβyⁱ].∫yds[1]-nls[:∫∂vβyⁱ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫∂vβyⁱ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫βγyⁱ].∫yds[1]-nls[:∫βγyⁱ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫βγyⁱ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			aerr = abs(nls2x[:∫γyⁱ].∫yds[1]-nls[:∫γyⁱ].∫yds[1]);
			rerr = aerr/abs(nls2x[:∫γyⁱ].∫yds[1]);
			if (aerr>prm[:atol][1]&&rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			#  solution cohort nodal errors
			myerrs!( (prm[:yˢrg][2]-prm[:yˢrg][1])*sol.yˢ.ys,
				 (prm[:yˢrg][2]-prm[:yˢrg][1])*sol2x.yˢ.ys,
				 yˢlow,yˢhigh;
				 rerr=yˢrerr,aerr=yˢaerr,
				 ngrp=yˢngrp,ntunit=yˢntunit );
			aerr = abs(sol.yˢ.ys[1]-sol2x.yˢ.ys[1]);
			rerr = aerr/abs(sol2x.yˢ.ys[1]);
			aerr *= (prm[:yˢrg][2]-prm[:yˢrg][1]);
			if (!myerrtst(yˢaerr,yˢrerr;prm=prm))||( aerr>prm[:atol][1]&&rerr>prm[:rtol][1] )
				Δt*=0.5;
				continue
			end

			myerrs!( (prm[:yᵛrg][2]-prm[:yᵛrg][1])*sol.yᵛ.ys,
				 (prm[:yᵛrg][2]-prm[:yᵛrg][1])*sol2x.yᵛ.ys,
				 yᵛlow,yᵛhigh;
				 rerr=yᵛrerr,aerr=yᵛaerr,
				 ngrp=yᵛngrp,ntunit=yᵛntunit );
			aerr = abs(sol.yᵛ.ys[1]-sol2x.yᵛ.ys[1]);
			rerr = aerr/abs(sol2x.yᵛ.ys[1]);
			aerr *= (prm[:yᵛrg][2]-prm[:yᵛrg][1]);
			if (!myerrtst(yᵛaerr,yᵛrerr;prm=prm))||( aerr>prm[:atol][1]&&rerr>prm[:rtol][1] )
				Δt*=0.5;
				continue
			end
			
			myerrs!( (prm[:yⁱrg][2]-prm[:yⁱrg][1])*sol.yⁱ.ys,
				 (prm[:yⁱrg][2]-prm[:yⁱrg][1])*sol2x.yⁱ.ys,
				 yⁱlow,yⁱhigh;
				 rerr=yⁱrerr,aerr=yⁱaerr,
				 ngrp=yⁱngrp,ntunit=yⁱntunit );
			aerr = abs(sol.yⁱ.ys[1]-sol2x.yⁱ.ys[1]);
			rerr = aerr/abs(sol2x.yⁱ.ys[1]);
			aerr *= (prm[:yⁱrg][2]-prm[:yⁱrg][1]);
			if (!myerrtst(yⁱaerr,yⁱrerr;prm=prm))||( aerr>prm[:atol][1]&&rerr>prm[:rtol][1] )
				Δt*=0.5;
				continue
			end

			aerr = abs(yʳ[1]-yʳ2x[1]);
			rerr = aerr/abs(yʳ2x[1]);
			if (aerr>prm[:atol][1])&&(rerr>prm[:rtol][1])
				Δt*=0.5;
				continue
			end

			flagfd = true;
		end

		# If any values were negative, set to 0. 
		@inbounds for i=1:nnd
			sol2x.yˢ.ys[i] = sol2x.yˢ.ys[i]>=0 ? sol2x.yˢ.ys[i] : 0.0;
			sol2x.yᵛ.ys[i] = sol2x.yᵛ.ys[i]>=0 ? sol2x.yᵛ.ys[i] : 0.0;
			sol2x.yⁱ.ys[i] = sol2x.yⁱ.ys[i]>=0 ? sol2x.yⁱ.ys[i] : 0.0;
		end

		# solution was accepted so prepare for next iteration
		#  store solution if at a downsample
		val = sol2x.t₀[1]-prm[:Trg][1];
		while (val >= prg+δprg)&&(pos<=nsmp)
			ysol[pos] = myinterp(taxis[pos],sol0,sol2x) |> (x->srefine(nndsmp,x));
			yʳsol[pos] = myinterp([sol0.t₀[1],sol2x.t₀[1]],[yʳ[1],yʳ2x[1]],taxis[pos]);
			pos += 1;
			prg += δprg;
			if flagprg
				printprg = prg/(prm[:Trg][2]-prm[:Trg][1]);
				println("Simulation progress: $printprg/1"); 
			end
		end
		sol0 = deepcopy(sol2x); nonlocals!(sol0;temp=nls,prm=prm); yʳ0[:] = deepcopy(yʳ2x);
		Δt *= 2;

		if ( (sol0.yˢ.∫yds[1]<-0.1)||(sol0.yᵛ.∫yds[1]<-0.1)||(sol0.yⁱ.∫yds[1]<-0.1)||
		      (sol0.yˢ.∫yds[1]+sol0.yᵛ.∫yds[1]+sol0.yⁱ.∫yds[1]+yʳ2x[1] > 1.075+prm[:ρ][1])||
		      (sol0.yˢ.∫yds[1]+sol0.yᵛ.∫yds[1]+sol0.yⁱ.∫yds[1]+yʳ2x[1] < 0.925+prm[:ρ][1]) ) 
			@warn "simulation failed at these tolerances"
			break
		end
	end

	prm[:nΔtfail][1] = nΔtfail;
	println("Number of times Δt was too small for refinement: $nΔtfail");
	return ysol,yʳsol
end
