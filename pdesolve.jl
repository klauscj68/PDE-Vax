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
		βy.ylvl[i] = β(ylvl.tlvl.nds[:,i],prm;case=:χτ)*
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
		λy.ylvl[i] = λ(ylvl.tlvl.nds[:,i],prm;case=:χτ)*
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
		Imαy.ylvl[i] = (1-α(ylvl.tlvl.nds[:,i],prm;case=:χτ))*
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
optional dictionary input for writing solutions and βy,λy,(1-α)y values

Note: Keep the EYSOL ∫yds fields empty and then euler! will internally 
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

	# Update the non-∂ nodal solution values for this Euler step 
	for i=2:nnd
		δτ = EYSOL[:yˢ].tlvl.nds[2,i] - YSOL[:yˢ].tlvl.nds[2,i];
		EYSOL[:yˢ].ys[i] = YSOL[:yˢ].ys[i-1] + δτ*
		                     -1/√(2)*YSOL[:yˢ].ys[i-1]*( YSOL[:λ].ys[i-1] + ∫βyⁱds );
		EYSOL[:yᵛ].ys[i] = YSOL[:yᵛ].ys[i-1] + δτ*
		                     -1/√(2)*YSOL[:yᵛ].ys[i-1]*( YSOL[:α].ys[i-1] + (1-YSOL[:α].ys[i-1])*∫βyⁱds );
		EYSOL[:yⁱ].ys[i] = YSOL[:yⁱ].ys[i-1] + δτ*
		                     -1/√(2)*YSOL[:γ].ys[i-1]*YSOL[:yⁱ].ys[i-1];
	end
	
	# Update the ∂-nodal solution value at (-t₀,0)
	EYSOL[:yˢ].ys[1] = 0.;
	EYSOL[:yᵛ].ys[1] = ∫λyˢds;
	EYSOL[:yⁱ].ys[1] = (∫yˢds + ∫Imαyᵛds)*∫βyⁱds;
	
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

# vaxsolver
"""
Integrate the vaccination system by an explicit Euler scheme with adaptive timestep
"""
function vaxsolver(prm::Dict{Symbol,Float64})
	# Define initial data
	∂YSOL = Dict{Symbol,Yℓvℓ}();

	
end
