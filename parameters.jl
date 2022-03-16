using SpecialFunctions
"""
Output the default model parameters as a Dict{Symbol,Vector{Float64}}
"""
function data()
	prm = Dict{Symbol,Vector{Float64}}();

	# domain
	#  Ranges should be the support of your initial data
	#  data! will adjust so that domains correctly contain
	#  support of solution over simulation time span. Needed
	#  to correctly compute the bdflow equations
	prm[:yˢrg₀] = [0.0,100.0*365];
	prm[:yᵛrg₀] = [0.0,0.0];
	prm[:yⁱrg₀] = [0.0,14.0];

	#  These are overwritten in data! to match actual simulation
	#  parameters
	prm[:yˢrg] = [NaN,NaN];
	prm[:yᵛrg] = [NaN,NaN];
	prm[:yⁱrg] = [NaN,NaN];

	prm[:Trg] = [0.0,138.0];

	# α parameters
	prm[:αL] = [14.0];
	prm[:αeff] = [0.92];

	# β parameters
	#  mean of 3.1
	prm[:βθ]=[3.5];
	prm[:βα]=[2.0];

	# γ parameters
	#  mean of 7.1
	prm[:γθ]=[8.0];
	prm[:γα]=[2.5];

	# λ parameters
	prm[:λ]=[1.5];

	# mass of t₀ infected
	prm[:ρ] = [0.0051*1e2];

	# spatial discretization
	prm[:nnd] = [500.0];

	# ode discretization
	prm[:atol]=[1e-6];
	prm[:rtol]=[1e-3];

	# parameters for how often and the res by which sol is stored
	prm[:dwnsmp]=[1.0];
	prm[:nndsmp]=[500.0];

	# normalization constants for fˢ,fⁱ. data! will mutate to correct
	prm[:fˢη]=[1.0];
	prm[:fⁱη]=[1.0];

	# placeholder for writing the ℓerr in abc
	prm[:ℓerr] = [NaN];

	return prm
end

"""
Adjust the prm data set for variables that depend on others. 
"""
function data!(prm::DSymVFl)
	nnd = prm[:nnd][1]|>ceil|>Int64;

	# Adjust the yrg's to correctly contain the support of initial data
	# during simulation time span
	prm[:yˢrg][1] = prm[:yˢrg₀][1];
	prm[:yᵛrg][1] = prm[:yᵛrg₀][1];
	prm[:yⁱrg][1] = prm[:yⁱrg₀][1];

	prm[:yˢrg][2] = prm[:yˢrg₀][2] + prm[:Trg][2];
	prm[:yᵛrg][2] = 31.0; #prm[:yᵛrg₀][2] + prm[:Trg][2];
	prm[:yⁱrg][2] = 31.0; #prm[:yⁱrg₀][2] + prm[:Trg][2];

	# Compute the normalization constants for initial conditions
	#  First reset to 1 so get correct new factor
	prm[:fˢη][1] = 1.0;
	prm[:fⁱη][1] = 1.0;

	#  Compute fˢ discretization and normalization
	tlvl = Tℓvℓ(0.0,nnd,prm[:yˢrg]);
	ys   = [fˢ(s;prm=prm) for s in tlvl.snds];
	ylvl = Yℓvℓ(tlvl,ys);
	∫line!(ylvl);
	prm[:fˢη][1] = 1/ylvl.∫yds[1];

	#  Compute fⁱ discretization and normalization
	tlvl = Tℓvℓ(0.0,nnd,prm[:yⁱrg]);
	ys   = [fⁱ(s;prm=prm) for s in tlvl.snds];
	ylvl = Yℓvℓ(tlvl,ys);
	∫line!(ylvl);
	prm[:fⁱη][1] = prm[:ρ][1]/ylvl.∫yds[1];

	# Compute the ρfⁱ(0) value
	val = fⁱ(0.0;prm=prm);

	# Compute the ∫βyⁱ at t=0
	tlvl = Tℓvℓ(0.0,nnd,prm[:yⁱrg]);
	ys = [β(s,0.0;prm=prm)*fⁱ(s;prm=prm) for s in tlvl.snds];
	ylvl = Yℓvℓ(tlvl,ys);
	∫line!(ylvl);
	∫val = ylvl.∫yds[1];

	ρ = val/∫val;

	# Compute the βθ value which matches and set dictionary to that
	# 1/(ηa^b) = ρ/a^b => 1/η^b=ρ => η=ρ^-1/b
	prm[:βθ][1] =ρ^(-1/prm[:βα][1])*prm[:βθ][1];
	
end

#%% Equation terms
"""
Evaluate the α equation terms for given choice of parameters
Note: ∂v = ∂s+∂t
"""
function α(s::Float64,t::Float64;prm::DSymVFl=data())
	val = s<=0.0 ? 0.0 : prm[:αeff][1]*exp(-prm[:αL][1]/14/s);

	return val
end
function ∂vα(s::Float64,t::Float64;prm::DSymVFl=data())
	val = s<=0.0 ? 0.0 : prm[:αeff][1]*exp(-prm[:αL][1]/14/s)*prm[:αL][1]/14/s^2;
end

function β(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:βα][1]/prm[:βθ][1]*(s/prm[:βθ][1])^(prm[:βα][1]-1);

	return val
end
function ∂vβ(s::Float64,t::Float64;prm::DSymVFl=data())
	return prm[:βα][1]/prm[:βθ][1]*(s/prm[:βθ][1])^(prm[:βα][1]-2.)*(prm[:βα][1]-1.)*1/prm[:βθ][1]

end

function γ(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:γα][1]/prm[:γθ][1]*(s/prm[:γθ][1])^(prm[:γα][1]-1);

	return val
end
function ∂vγ(s::Float64,t::Float64;prm::DSymVFl=data())
	return prm[:γα][1]/prm[:γθ][1]*(s/prm[:γθ][1])^(prm[:γα][1]-2.)*(prm[:γα][1]-1.)*1/prm[:γθ][1]
end
function λ(s::Float64,t::Float64;prm::DSymVFl=data())
	# Formula: b₁(s)*{ b₂(s)*[(1-σ)*λ₁(t)+(σ-1)*λ₂(t)] + λ₂(t) }
	syr = s/365;
	# Two week rollout before steady state hazard
	if (syr<=20.0) # Age 20yrs or before Dec 15th (starting Oct 1)
		return 0.0
	elseif (syr>20.0)&&(syr<70.0)
		b₁ = exp(-2/(syr-20));
		b₂ = exp(-2/(70-syr))/0.83;

		λ₁ = 0.002*mynrm(105.0,10.0,t) + 0.002*mynrm(140.0,7.5,t);
		λ₂ = 0.001*mynrm(93.0,6.0,t) + 0.013*mynrm(130.0,12.5,t);
		σ = (70.0-syr)/50;

		return b₁*( b₂*((1-σ)*λ₁+(σ-1)*λ₂) + λ₂ )
	else
		b₁ = exp(-2/(syr-20));
		λ₂ = 0.001*mynrm(93.0,6.0,t) + 0.013*mynrm(130.0,12.5,t);

		return b₁*λ₂
	end
end
function ∂vλ(s::Float64,t::Float64;prm::DSymVFl=data())
	syr = s/365;
	if (syr<=20.0)
		return 0.0
	elseif (syr>20.0)&&(syr<70.0)
		b₁ = exp(-2/(syr-20));
		b₂ = exp(-2/(70-syr))/0.83;

		λ₁ = 0.002*mynrm(105.0,10.0,t) + 0.002*mynrm(140.0,7.5,t);
		λ₂ = 0.001*mynrm(93.0,6.0,t) + 0.013*mynrm(130.0,12.5,t);
		σ = (70.0-syr)/50;

		∂b₁ = exp(-2/(syr-20))*(2/(syr-20)^2)/365;
		∂b₂ = exp(-2/(70-syr))/0.83*(-2/(70-syr)^2)/365;

		∂λ₁ = 0.002*∂mynrm(105.0,10.0,t) + 0.002*∂mynrm(140.0,7.5,t);
		∂λ₂ = 0.001*∂mynrm(93.0,6.0,t) + 0.013*∂mynrm(130.0,12.5,t);
		∂σ = -1/50/365;

		val = ( ∂b₁*(b₂*((1-σ)*λ₁+(σ-1)*λ₂) + λ₂) + 
		       b₁*( ∂b₂*((1-σ)*λ₁+(σ-1)*λ₂)
			   + b₂*(-∂σ*λ₁ + (1-σ)*∂λ₁ + ∂σ*λ₂ + (σ-1)*∂λ₂) + ∂λ₂
 		          )
		      );

		return val
	else
		b₁ = exp(-2/(syr-20));
		λ₂ = 0.001*mynrm(93.0,6.0,t) + 0.013*mynrm(130.0,12.5,t);

		∂b₁ = exp(-2/(syr-20))*(2/(syr-20)^2)/365;
		∂λ₂ = 0.001*∂mynrm(93.0,6.0,t) + 0.013*∂mynrm(130.0,12.5,t);

		val = ∂b₁*λ₂ + b₁*∂λ₂
		return val
	end
end

#%% Initial data
function fˢ(s::Float64;prm::DSymVFl=data())
	syr = s/365;
	if (s<=0.0)||(syr>=100)
		return 0.0
	else
		return prm[:fˢη][1]*exp(-1/s)*exp(-6/(100-s))
	end
end
function fᵛ(s::Float64;prm::DSymVFl=data())
	return 0.0
end
function fⁱ(s::Float64;prm::DSymVFl=data())
	if (s>=14.0)
		return 0.0
	else
		return prm[:ρ][1]*prm[:fⁱη][1]*exp(-1e-3/(14-s))
	end
end
