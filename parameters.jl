## File storing parameters and model equation terms
# Should also have loaded auxilliary.jl

#%% Model scalar parameters
# data
"""
Dictionary storing model scalar parameters
"""
function data()
	prm = Dict{Symbol,Float64}();
	
	# Domain
	#  Largest age
	prm[:L] = 10.;

	#  Largest time
	prm[:T] = 10.;	

	# Epidemic
	#  Initial fraction vaccinated
	prm[:ρ] = .3;

	# Hazard rates
	#  β
	prm[:βa] = 1.;
	prm[:βb] = 1.;
	prm[:βη] = 1.; # mutated by ∂YSOL! to ensure continuity at (0,0)

	#  γ
	prm[:γa] = 1.;
	prm[:γb] = 1.;

	# Numerical discretization
	#  Number nodes within each [t=t₀] set
	prm[:nnd] = 3.;

	#  Final t-time of integration
	prm[:tfin] = 3.;

	#   t-downsample used to save the solution along integration
	prm[:δt] = 1.;

	#   tolerances for ode integration and maximum permitted 
	#   integration step
	prm[:atol] = 1e-6;
	prm[:rtol] = 1e-3;
	prm[:δtmax] = prm[:δt];

	return prm
	
end

#%% Model function parameters
# λ
"""
Evaluate λ 
Defaults to λ(s,t) but optional case argument can be used to instead
compute λ(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function λ(pt::VecVw,prm::Dict{Symbol,Float64};case::Symbol=:st)
	
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of λ(s,t) given here
		#  Note: need λ(s,0)≡0 for BC's
		val = (t<=10. ? t : 10.);
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = λ(newpt,prm;case=:st);
	else
		error("not valid λ-eval case");
	end

	return val
end

# β
"""
Evaluate β
Defaults to β(s,t) but optional case argument can be used to instead
compute β(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function β(pt::VecVw,prm::Dict{Symbol,Float64};case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of β(s,t) given here	
		val = prm[:βb]/prm[:βa]*real( (s/prm[:βa]+0im)^(prm[:βb]-1.) );
		val *= prm[:βη]; # used by ∂YSOL! to enforce BC 
		                 # continuity at (0,0)
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = β(newpt,prm;case=:st);
	else
		error("not valid β-eval case");
	end

	return val
end

# α
"""
Evaluate α
Defaults to α(s,t) but optional case argument can be used to instead
compute α(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function α(pt::VecVw,prm::Dict{Symbol,Float64};case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of α(s,t) given here	
		val = .92/14*(s<=14 ? s : 14);
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = α(newpt,prm;case=:st);
	else
		error("not valid α-eval case");
	end

	return val
end

# γ
"""
Evaluate γ
Defaults to γ(s,t) but optional case argument can be used to instead
compute γ(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function γ(pt::VecVw,prm::Dict{Symbol,Float64};case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of γ(s,t) given here	
		val = (s != 0 ? 
		        prm[:γb]/prm[:γa]*real( (s/prm[:γa]+0im)^(prm[:γb]-1.) )
			: 0.);
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = γ(newpt,prm;case=:st);
	else
		error("not valid γ-eval case");
	end

	return val
end

# fˢ
"""
Evaluate fˢ
Defaults to fˢ(s,t) but optional case argument can be used to instead
compute fˢ(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function fˢ(pt::VecVw,prm::Dict{Symbol,Float64};case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of fˢ given here	
		val = 1. /prm[:L];
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = fˢ(newpt,prm;case=:st);
	else
		error("not valid fˢ-eval case");
	end

	return val
end

#fⁱ
"""
Evaluate fⁱ
Defaults to fⁱ(s,t) but optional case argument can be used to instead
compute fⁱ(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function fⁱ(pt::VecVw,prm::Dict{Symbol,Float64};case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of fⁱ given here	
		val = 1. /prm[:L];
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = fⁱ(newpt,prm;case=:st);
	else
		error("not valid fⁱ-eval case");
	end

	return val
end

#%% Auxilliary methods for sampling function parameters across batch of sample points
function λ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = λ(pt,prm;case=case);
	end

	return val
end
function β(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = β(pt,prm;case=case);
	end

	return val
end
function α(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = α(pt,prm;case=case);
	end

	return val
end
function γ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = γ(pt,prm;case=case);
	end

	return val
end
function fˢ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = fˢ(pt,prm;case=case);
	end

	return val
end
function fⁱ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = fⁱ(pt,prm;case=case);
	end

	return val
end
