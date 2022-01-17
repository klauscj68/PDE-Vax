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

	prm[:Trg] = [0.0,0.5];

	# α parameters
	prm[:αL] = [14.0];
	prm[:αeff] = [0.92];

	# β parameters
	prm[:βa]=[1.0];
	prm[:βb]=[3.0];

	# γ parameters
	prm[:γa]=[1.0];
	prm[:γb]=[2.0];

	# λ parameters
	prm[:λ]=[0.5];

	# mass of t₀ infected
	prm[:ρ] = [0.3];

	# spatial discretization
	prm[:nnd] = [50.0];

	# ode discretization
	prm[:atol]=[1e-6];
	prm[:rtol]=[1e-3];

	# parameters for how often and the res by which sol is stored
	prm[:dwnsmp]=[0.1];
	prm[:nndsmp]=[5.0];

	return prm
end

#%% Equation terms
"""
Evaluate the α equation terms for given choice of parameters
Note: ∂v = ∂s+∂t
"""
function α(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:αeff][1]*(s<prm[:αL][1] ? s : prm[:αL][1]);

	return val
end
function ∂vα(s::Float64,t::Float64;prm::DSymVFl=data())
	if s>prm[:αL][1]
		return 0.0
	else
		return prm[:αeff][1]/prm[:αL][1]
	end
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
	return prm[:λ][1]
end
function ∂vλ(s::Float64,t::Float64;prm::DSymVFl=data())
	return 0.0
end

#%% Initial data
function fˢ(s::Float64;prm::DSymVFl=data())
	return 1/(100*365);
end
function fᵛ(s::Float64;prm::DSymVFl=data())
	return 0.0
end
function fⁱ(s::Float64;prm::DSymVFl=data())
	return prm[:ρ][1]*1/14;
end
