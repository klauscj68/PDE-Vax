# Routines used to integrate the vaccination pde system
#%% Ancillary routines for computing system integrands
# βy!
"""
Compute βy at the nodes along Tℓvℓ. Routine mutates the optional
argument βy.
"""
function βy!(ylvl::Yℓvℓ,prm::Dict{Symbol,Float64};
	     βy::Yℓvℓ=Yℓvℓ(undef))

	if isnan(βy.ys[1])
		βy = deepcopy(ylvl);
		flagrt = true;
	else
		flagrt = false;
	end
	
	for i=1:ylvl.tlvl.nnd
		βy.ys[i] = β(ylvl.tlvl.nds[:,i],prm;case=:χτ)*
		             ylvl.ys[i];
	end

	if flagrt
		return βy
	end
end

# λy!
"""
Compute λy at the nodes along Tℓvℓ. Routine mutates the optional 
argument λy
"""
function λy!(ylvl::Yℓvℓ,prm::Dict{Symbol,Float64};
	     λy::Yℓvℓ=Yℓvℓ(undef))

	if isnan(λy.ys[1])
		λy = deepcopy(ylvl);
		flagrt = true;
	else
		flagrt = false;
	end

	for i=1:ylvl.tlvl.nnd
		λy.ys[i] = λ(ylvl.tlvl.nds[:,i],prm;case=:χτ)*
		             ylvl.ys[i];
	end

	if flagrt
		return λy
	end
end
		

# Imαy!
"""
Compute (1-α)y at the nodes along Tℓvℓ. Routine mutates the optional
argument Imαy
"""
function Imαy!(ylvl::Yℓvℓ,prm::Dict{Symbol,Float64};
	       Imαy::Yℓvℓ=Yℓvℓ(undef))
	
	if isnan(Imαy.ys[1])
		Imαy = deepcopy(ylvl);
		flagrt = true;
	else
		flagrt = false;
	end

	for i=1:ylvl.tlvl.nnd
		Imαy.ys[i] = (1-α(ylvl.tlvl.nds[:,i],prm;case=:χτ))*
			      ylvl.ys[i];
	end

	if flagrt
		return Imαy
	end
end

#%% System integrator
# euler!
"""
Compute an explicit Euler step of the PDE vaccination system, mutates the
optional dictionary input for writing solutions and λ,α,γ and 
βy,λy,(1-α)y values

Note: Keep the EYSOL ∫yds fields NaN and then euler! will internally 
      mutate those integrals to their correct values when in later iteration
      it is passed as YSOL
"""
function euler!(δt::Float64,YSOL::Dict{Symbol,Yℓvℓ},prm::Dict{Symbol,Float64};
		EYSOL::Dict{Symbol,Yℓvℓ}=Dict{Symbol,Yℓvℓ}())
	
	nnd = YSOL[:yˢ].tlvl.nnd;
	t₀ = YSOL[:yˢ].tlvl.t₀[1];	
	if isempty(EYSOL)
		EYSOL = deepcopy(YSOL);
		flagrt = true;
	else
		flagrt = false;
	end
	
	# Compute needed integrals for forward Euler step
	#  Note the integral value is constant within [t=t₀]
	if isnan(YSOL[:βyⁱ].∫yds[1])
		∫βyⁱds = ∫line(YSOL[:βyⁱ]);
		YSOL[:βyⁱ].∫yds[1] = ∫βyⁱds;
	else
		∫βyⁱds = YSOL[:βyⁱ].∫yds[1];
	end

	if isnan(YSOL[:λyˢ].∫yds[1])
		∫λyˢds = ∫line(YSOL[:λyˢ]);
		YSOL[:λyˢ].∫yds[1] = ∫λyˢds;
	else
		∫λyˢds = YSOL[:λyˢ].∫yds[1];
	end

	if isnan(YSOL[:yˢ].∫yds[1])
		∫yˢds = ∫line(YSOL[:yˢ]);
		YSOL[:yˢ].∫yds[1] = ∫yˢds;
	else
		∫yˢds = YSOL[:yˢ].∫yds[1];
	end

	if isnan(YSOL[:Imαyᵛ].∫yds[1])
		∫Imαyᵛds = ∫line(YSOL[:Imαyᵛ]);
		YSOL[:Imαyᵛ].∫yds[1] = ∫Imαyᵛds;
	else
		∫Imαyᵛds = YSOL[:Imαyᵛ].∫yds[1];
	end

	# Advance to new tlvl's for this Euler step
	Tℓvℓ!(t₀+δt,YSOL[:yˢ].tlvl,EYSOL[:yˢ].tlvl); Tℓvℓ!(t₀+δt,YSOL[:λ].tlvl,EYSOL[:λ].tlvl);
	Tℓvℓ!(t₀+δt,YSOL[:yᵛ].tlvl,EYSOL[:yᵛ].tlvl); Tℓvℓ!(t₀+δt,YSOL[:α].tlvl,EYSOL[:α].tlvl);
	Tℓvℓ!(t₀+δt,YSOL[:yⁱ].tlvl,EYSOL[:yⁱ].tlvl); Tℓvℓ!(t₀+δt,YSOL[:γ].tlvl,EYSOL[:γ].tlvl);
	
	Tℓvℓ!(t₀+δt,YSOL[:βyⁱ].tlvl,EYSOL[:βyⁱ].tlvl);
	Tℓvℓ!(t₀+δt,YSOL[:λyˢ].tlvl,EYSOL[:λyˢ].tlvl);
	Tℓvℓ!(t₀+δt,YSOL[:Imαyᵛ].tlvl,EYSOL[:Imαyᵛ].tlvl);

	# Update the ∂-nodal solution value at (χ,τ) = (-t₀,0)
	EYSOL[:yˢ].ys[1] = 0.;
	EYSOL[:yᵛ].ys[1] = ∫λyˢds;
	EYSOL[:yⁱ].ys[1] = (∫yˢds + ∫Imαyᵛds)*∫βyⁱds;
	
	# Update the non-∂ nodal solution values for this Euler step 
	χsY = @view YSOL[:yˢ].tlvl.nds[1,:];
	@inbounds for i=2:nnd
		# Find value of solution directly underneath 
		nd = @view EYSOL[:yˢ].tlvl.nds[:,i];
		if nd[1] >= χsY[1]
			# Node is overtop the prior t level
			nd0 = γℓvℓ(t₀,nd[1]);
			δτ = nd[2] - nd0[2];

			yˢ₀ = myinterp(χsY,YSOL[:yˢ].ys,nd[1]);
			yᵛ₀ = myinterp(χsY,YSOL[:yᵛ].ys,nd[1]);
			yⁱ₀ = myinterp(χsY,YSOL[:yⁱ].ys,nd[1]);

		else
			# Node is overtop the ∂
			nd0 = [nd[1],0.];
			δτ = EYSOL[:yˢ].tlvl.nds[2,1];
			
			yˢ₀ = 0.;
			yᵛ₀ = myinterp([ EYSOL[:yᵛ].tlvl.nds[1,1],χsY[1] ],
				       [ ∫λyˢds,YSOL[:yᵛ].ys[1] ],
				       nd[1]);
			yⁱ₀ = myinterp([ EYSOL[:yⁱ].tlvl.nds[1,1],χsY[1] ],
				       [ (∫yˢds + ∫Imαyᵛds)*∫βyⁱds,YSOL[:yⁱ].ys[1] ],
				       nd[1]);

		end

		λ₀ = λ(nd0,prm;case=:χτ);
		α₀ = α(nd0,prm;case=:χτ);
		γ₀ = γ(nd0,prm;case=:χτ);

		EYSOL[:yˢ].ys[i] = yˢ₀ - 1/√(2)*yˢ₀*( λ₀ + ∫βyⁱds )*δτ;
		EYSOL[:yᵛ].ys[i] = yᵛ₀ - 1/√(2)*yᵛ₀*( α₀ + (1-α₀)*∫βyⁱds )*δτ;
		EYSOL[:yⁱ].ys[i] = yⁱ₀ - 1/√(2)*γ₀*yⁱ₀*δτ;
	end	
	
	# Update the intermediate quantities
	βy!(EYSOL[:yⁱ],prm; βy=EYSOL[:βyⁱ]);
	λy!(EYSOL[:yˢ],prm; λy=EYSOL[:λyˢ]);
	Imαy!(EYSOL[:yᵛ],prm; Imαy=EYSOL[:Imαyᵛ]);

	for i=1:nnd
		EYSOL[:λ].ys[i] = λ(EYSOL.tlvl.nds[i],prm;case=:χτ);
		EYSOL[:α].ys[i] = α(EYSOL.tlvl.nds[i],prm;case=:χτ);
		EYSOL[:γ].ys[i] = γ(EYSOL.tlvl.nds[i],prm;case=:χτ);
	end

	if flagrt
		return EYSOL
	end

end

# ∂YSOL
"""
Routine to define the boundary data at the [t=0.] level
"""
function ∂YSOL(prm::Dict{Symbol,Float64})
	# Define initial geometry
	tlvl = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:L],Int64(prm[:nnd]))) );
	
	# Define [t=0.] ∂-values
	YSOL = Dict{Symbol,Yℓvℓ}();

	yˢ = Vector{Float64}(undef,prm[:nnd]);
	yᵛ = zeros(prm[:nnd]);
	yⁱ = Vector{Float64}(undef,prm[:nnd]);
	λs = Vector{Float64}(undef,prm[:nnd]);
	αs = Vector{Float64}(undef,prm[:nnd]);
	γs = Vector{Float64}(undef,prm[:nnd]);
	for i=1:prm[:nnd]
		yˢ[i] = fˢ(@view tlvl.nds[:,i],prm;case=:χτ);
		yⁱ[i] = prm[:ρ]*fⁱ(@view tlvl.nds[:,i],prm;case=:χτ);

		λs[i] = λ(@view tlvl.nds[:,i],prm;case=:χτ);
		αs[i] = α(@view tlvl.nds[:,i],prm;case=:χτ);
		γs[i] = γ(@view tlvl.nds[:,i],prm;case=:χτ);
	end
	YSOL[:yˢ] = Yℓvℓ(tlvl,yˢ);
	YSOL[:yᵛ] = Yℓvℓ(tlvl,yᵛ);
	YSOL[:yⁱ] = Yℓvℓ(tlvl,yⁱ);
	
	#  βyⁱ,λyˢ,Imαyᵛ
	YSOL[:βyⁱ] = βy!(YSOL[:yⁱ],prm);
	YSOL[:λyˢ] = λy!(YSOL[:yˢ],prm);
	YSOL[:Imαyᵛ] = Imαy!(YSOL[:yᵛ],prm);

	return YSOL
end

# vaxsolver
"""
Integrate the vaccination system by an explicit Euler scheme with adaptive timestep
"""
function vaxsolver(prm::Dict{Symbol,Float64})
	# Intialize values and vectors for storing solution 
	taxis = convert(Vector,0.: prm[:δt] : prm[:tfin]);
	ntdwn = length(taxis);
	SOL = Vector{Dict{Symbol,Yℓvℓ}}(undef,ntdwn);
	
	# Setup ∂-[t=0] data
	SOL[1] = ∂YSOL(prm);

	# Adaptive Euler step requires 4 Yℓvℓ mem locs for writing output
	ynow = deepcopy(SOL[1]); ynext = deepcopy(ynow);
	ymid = deepcopy(ynow); y2xmid = deepcopy(ynow);
	
	pos = 2; # indicates which taxis value is next to surpass
	δt = prm[:δt]; # initial guess of adaptive Euler step
	yaerr = [0.,0.,0.]; # stores the ODE solver absolute error
	yrerr = [0.,0.,0.]; # stores the ODE solver relative error
	while pos <= ntdwn
		# Compute the full step
		euler!(δt,ynow,prm;EYSOL=ynext);

		# Compute the two half-steps
		euler!(.5*δt,ynow,prm;EYSOL=ymid);
		euler!(.5*δt,ymid,prm;EYSOL=y2xmid);

		# Compute the abs errors
		yaerr[1] = maximum(abs.(y2xmid[:yˢ].ys-ynext[:yˢ].ys));
		yaerr[2] = maximum(abs.(y2xmid[:yᵛ].ys-ynext[:yᵛ].ys));
		yaerr[3] = maximum(abs.(y2xmid[:yⁱ].ys-ynext[:yⁱ].ys));

		# Compute the rel errors
		yrerr[1] = yaerr[1]/(sum(abs.(ynow[:yˢ]))/prm[:nnd]);
		yrerr[2] = yaerr[2]/(sum(abs.(ynow[:yᵛ]))/prm[:nnd]);
		yrerr[3] = yaerr[3]/(sum(abs.(ynow[:yⁱ]))/prm[:nnd]);


		# Act according to accepting or addapting the t-step
		flagδt = (yaerr .<= prm[:atol]).|(yrerr .<= prm[:rtol]);
		flag = flagδt[1]&&flagδt[2]&&flagδt[3];

		if !flag
			# Error tolerances not attained so adapt and iterate
			δt *= .5;
			continue
		end
		
		# Error tolerances attained so accept	
		tnow = ynow.tlvl.t₀[1]; tnext = y2xmid.tlvl.t₀[1];
		if tnext >= taxis[pos]	
			# We've passed a downsample t-value and will now store
			posnext = (tnext < taxis[end] ? 
				    myfindfirst(taxis,tnext) : ntdwn + 1 );

			# Interpolate the inbetween values
			for i=pos:posnext-1
				SOL[i][:yˢ] = myinterp([tnow,tnext],[ynow[:yˢ],y2xmid[:yˢ]],taxis[i]);
				SOL[i][:yᵛ] = myinterp([tnow,tnext],[ynow[:yᵛ],y2xmid[:yᵛ]],taxis[i]);
				SOL[i][:yⁱ] = myinterp([tnow,tnext],[ynow[:yⁱ],y2xmid[:yⁱ]],taxis[i]);

				SOL[i][:λ] = myinterp([tnow,tnext],[ynow[:λ],y2xmid[:λ]],taxis[i]);
				SOL[i][:α] = myinterp([tnow,tnext],[ynow[:α],y2xmid[:α]],taxis[i]);
				SOL[i][:γ] = myinterp([tnow,tnext],[ynow[:γ],y2xmid[:γ]],taxis[i]);

				SOL[i][:βyⁱ] = myinterp([tnow,tnext],[ynow[:βyⁱ],y2xmid[:βyⁱ]],taxis[i]);
				SOL[i][:λyˢ] = myinterp([tnow,tnext],[ynow[:λyˢ],y2xmid[:λyˢ]],taxis[i]);
				SOL[i][:Imαyᵛ] = myinterp([tnow,tnext],[ynow[:Imαyᵛ],y2xmid[:Imαyᵛ]],taxis[i]);
			end

			# Define new position
			pos = posnext;
	
		end
		
		for key in keys(ynow)
			ynow[key].tlvl.t₀[:] = y2xmid[key].tlvl.t₀;
			ynow[key].tlvl.nds[:,:] = y2xmid[key].tlvl.nds;
			ynow[key].tlvl.χrg[:] = y2xmid[key].tlvl.χrg;
			ynow[key].tlvl.τrg[:] = y2xmid[key].tlvl.τrg;
			ynow[key].ys[:] = y2xmid[key].ys;
			ynow[key].∫yds[:] = [NaN];
		end

		# Try a larger time step and continue iteration
		δt = minimum([2*δt,prm[:δtmax]]);
	end

	return taxis,SOL
end
