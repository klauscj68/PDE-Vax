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

	prm[:Trg] = [0.0,65.0];

	# α parameters
	prm[:αL] = [14.0];
	prm[:αeff] = [0.940723];

	# β parameters
	#  mean of 3.1
	prm[:βθ]=[13.5087];
	prm[:βα]=[18.6117];

	# γ parameters
	#  mean of 7.1
	prm[:γθ]=[14.5368];
	prm[:γα]=[7.96798];

	# λ parameters
	prm[:λ]=[1.5];

	# mass of t₀ infected
	prm[:ρ] = [0.0496982];

	# spatial discretization
	prm[:nnd] = [2500.0];

	# ode discretization
	prm[:atol]=[1e-5];
	prm[:rtol]=[1e-3];
	prm[:Δtmin]=[1e-5];

	# parameters for how often and the res by which sol is stored
	prm[:dwnsmp]=[1.0];
	prm[:nndsmp]=[2500.0];

	# normalization constants for fˢ,fⁱ. data! will mutate to correct
	prm[:fˢη]=[1.0];
	prm[:fⁱη]=[1.0];

	# placeholder for writing the ℓerr in abc
	prm[:ℓerr] = [NaN];

	# placeholder for storing how many times solver hit the Δtmin
	prm[:nΔtfail] = [NaN];
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
	prm[:yᵛrg][2] = 22.0; #prm[:yᵛrg₀][2] + prm[:Trg][2];
	prm[:yⁱrg][2] = 22.0; #prm[:yⁱrg₀][2] + prm[:Trg][2];

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
	
	return prm[:αeff][1]*Ηδρ(s-7.0;δ=7.0,ρ=1.0)
end
function ∂vα(s::Float64,t::Float64;prm::DSymVFl=data())
	
	return prm[:αeff][1]*∂Ηδρ(s-7.0;δ=7.0,ρ=1.0)
end

function β₀(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:βα][1]/prm[:βθ][1]*(s/prm[:βθ][1])^(prm[:βα][1]-1);

	return val
end
function β(s::Float64,t::Float64;prm::DSymVFl=data(),
		                 M::Float64=50.0)
	sₘ = prm[:βθ][1]*( prm[:βθ][1]/prm[:βα][1]*M )^( 1/(prm[:βα][1]-1) );
	val = β₀(s,t;prm=prm);

	val = val*Ηδρ(sₘ-s)+M*Ηδρ(s-sₘ);
	return val
end
function ∂vβ₀(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:βα][1]/prm[:βθ][1]*(s/prm[:βθ][1])^(prm[:βα][1]-2.)*(prm[:βα][1]-1.)*1/prm[:βθ][1];

	return val
end
function ∂vβ(s::Float64,t::Float64;prm::DSymVFl=data(),
	                           M::Float64=50.0)	
	sₘ = prm[:βθ][1]*( prm[:βθ][1]/prm[:βα][1]*M )^( 1/(prm[:βα][1]-1) );
	val = ∂vβ₀(s,t;prm=prm);

	val = val*Ηδρ(sₘ-s)-β₀(s,t;prm=prm)*∂Ηδρ(sₘ-s)+M*∂Ηδρ(s-sₘ)
	return val
end

function γ₀(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:γα][1]/prm[:γθ][1]*(s/prm[:γθ][1])^(prm[:γα][1]-1);

	return val
end
function γ(s::Float64,t::Float64;prm::DSymVFl=data(),
		                 M::Float64=50.0)
	sₘ = prm[:γθ][1]*( prm[:γθ][1]/prm[:γα][1]*M )^( 1/(prm[:γα][1]-1) );
	val = γ₀(s,t;prm=prm);

	val = val*Ηδρ(sₘ-s)+M*Ηδρ(s-sₘ);
	return val
end
function ∂vγ₀(s::Float64,t::Float64;prm::DSymVFl=data())
	val = prm[:γα][1]/prm[:γθ][1]*(s/prm[:γθ][1])^(prm[:γα][1]-2.)*(prm[:γα][1]-1.)*1/prm[:γθ][1];

	return val
end
function ∂vγ(s::Float64,t::Float64;prm::DSymVFl=data(),
	                           M::Float64=50.0)
	sₘ = prm[:γθ][1]*( prm[:γθ][1]/prm[:γα][1]*M )^( 1/(prm[:γα][1]-1) );
	val = ∂vγ₀(s,t;prm=prm);

	val = val*Ηδρ(sₘ-s)-γ₀(s,t;prm=prm)*∂Ηδρ(sₘ-s)+M*∂Ηδρ(s-sₘ)
	return val
end

function λ(s::Float64,t::Float64;prm::DSymVFl=data())
	syr = s/365;
	val = ζδρ(syr,20.0,59.0)*( 0.002*mynrm(t,60.0,10.0)+0.002*mynrm(t,95.0,7.5) );
	val += ζδρ(syr,60.0,69.0)*( 0.006*mynrm(t,85.0,12.5) );
	val += ζδρ(syr,70.0,95.0)*( 0.013*mynrm(t,85.0,12.5) );

	return val
end
function ∂vλ(s::Float64,t::Float64;prm::DSymVFl=data())
	syr = s/365; η=1/365;
	val = η*∂ζδρ(syr,20.0,59.0)*( 0.002*mynrm(t,60.0,10.0)+0.002*mynrm(t,95.0,7.5) );
	val += ζδρ(syr,20.0,59.0)*( 0.002*∂mynrm(t,60.0,10.0)+0.002*∂mynrm(t,95.0,7.5) );

	val += η*∂ζδρ(syr,60.0,69.0)*( 0.006*mynrm(t,85.0,12.5) );
	val += ζδρ(syr,60.0,69.0)*( 0.006*∂mynrm(t,85.0,12.5) );
	
	val += η*∂ζδρ(syr,70.0,95.0)*( 0.013*mynrm(t,85.0,12.5) );
	val += ζδρ(syr,70.0,95.0)*( 0.013*∂mynrm(t,85.0,12.5) );

	return val
end

#%% Initial data
function fˢ(s::Float64;prm::DSymVFl=data())
	syr = s/365;
	val = prm[:fˢη][1]*ζδρ(syr,2.0,98.0;δ=0.25,ρ=0.25);

	return val
end
function fᵛ(s::Float64;prm::DSymVFl=data())
	return 0.0
end
function fⁱ(s::Float64;prm::DSymVFl=data())
	val = prm[:ρ][1]*prm[:fⁱη][1]*Ηδρ(13.5-s;δ=0.25,ρ=0.1)
	
	return val
end
