using SpecialFunctions
"""
Output the default model parameters as a Dict{Symbol,Vector{Float64}}
"""
function data()
	prm = Dict{Symbol,Vector{Float64}}();

	# domain
	prm[:yˢrg] = [0.0,100.0*365];
	prm[:yᵛrg] = [0.0,2*31.0];
	prm[:yⁱrg] = [0.0,14.0];

	prm[:Trg] = [0.0,31.0];

	# α parameters
	prm[:αL] = [14.0];
	prm[:αeff] = [0.92];

	# β parameters
	#  mean of 3.1
	prm[:βa]=[3.5];
	prm[:βb]=[2.0];

	# γ parameters
	#  mean of 7.1
	prm[:γa]=[8.0];
	prm[:γb]=[2.5];

	# λ parameters
	prm[:λ]=[1.5];

	# mass of t₀ infected
	prm[:ρ] = [0.0051];

	# spatial discretization
	prm[:nnd] = [500.0];

	# ode discretization
	prm[:atol]=[1e-6];
	prm[:rtol]=[1e-4];

	# parameters for how often and the res by which sol is stored
	prm[:dwnsmp]=[0.1];
	prm[:nndsmp]=[100.0];

	return prm
end

"""
Adjust the prm data set for variables that depend on others
"""
function data!(prm::DSymVFl)
	# Compute the ρfⁱ(0) value
	val = fⁱ(0.0;prm=prm);

	# Compute the ∫βyⁱ at t=0
	nnd = prm[:nnd][1]|>ceil|>Int64;
	tlvl=Tℓvℓ(0.0,nnd,prm[:yⁱrg]);
	ys = [β(s,0.0;prm=prm)*fⁱ(s;prm=prm) for s in tlvl.snds];
	ylvl = Yℓvℓ(tlvl,ys);
	∫line!(ylvl);
	∫val = ylvl.∫yds[1];

	ρ = val/∫val;

	# Compute the βa value which matches and set dictionary to that
	# 1/(ηa^b) = ρ/a^b => 1/η^b=ρ => η=ρ^-1/b
	prm[:βa][1] =ρ^(-1/prm[:βb][1])*prm[:βa][1];
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
	val = prm[:βb][1]/prm[:βa][1]*(s/prm[:βa][1])^(prm[:βb][1]-1);

	return val
end
function ∂vβ(s::Float64,t::Float64;prm::DSymVFl=data())
	return prm[:βb][1]/prm[:βa][1]*(s/prm[:βa][1])^(prm[:βb][1]-2.)*(prm[:βb][1]-1.)*1/prm[:βa][1]

end

function γ(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:γb][1]/prm[:γa][1]*(s/prm[:γa][1])^(prm[:γb][1]-1);

	return val
end
function ∂vγ(s::Float64,t::Float64;prm::DSymVFl=data())
	return prm[:γb][1]/prm[:γa][1]*(s/prm[:γa][1])^(prm[:γb][1]-2.)*(prm[:γb][1]-1.)*1/prm[:γa][1]
end
function λ(s::Float64,t::Float64;prm::DSymVFl=data())
	# Two week rollout before steady state hazard
	if (s<=25550.0)||(t<=0.0)
		return 0.0
	else
		return prm[:λ][1]*exp(-1/t)*exp(-365/(s-70))
	end
end
function ∂vλ(s::Float64,t::Float64;prm::DSymVFl=data())
	if (s<=25550.0)||(t<=0.0)
		return 0.0
	else
		return prm[:λ][1]*( (1/t^2*exp(-1/t))*exp(-365/(s-70))
				     + exp(-1/t)*(365/(s-70)^2*exp(-365/(s-70))) )
	end
end

#%% Initial data
function fˢ(s::Float64;prm::DSymVFl=data())
	if s<=0.0
		return 0.0
	else
		return 2.884922631393703e-5*exp(-365/s)
	end
end
function fᵛ(s::Float64;prm::DSymVFl=data())
	return 0.0
end
function fⁱ(s::Float64;prm::DSymVFl=data())
	return prm[:ρ][1]*1/14;
end
