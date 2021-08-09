## File storing parameters and model equation terms
# Should also have loaded auxilliary.jl

#%% Model scalar parameters
# data
"""
Dictionary storing model scalar parameters
"""
function data()
	D = Dict{Symbol,Float64}();
	
	# Domain
	#  Largest age
	D[:L] = 85.;

	#  Largest time
	D[:T] = 85.;	

	# Epidemic
	#  Initial fraction vaccinated
	D[:ρ] = .3;

	# Numerical discretization
	#  Number nodes within each [t=t₀] set
	D[:nnd] = 3.;

	#  Final τ-time of integration
	D[:tfin] = .2;

	#   τ-downsample used to save the solution along integration
	D[:δt] = .1;

	#   tolerances for ode integration and maximum permitted 
	#   integration step
	D[:atol] = 1e-6;
	D[:rtol] = 1e-3;
	D[:δtmax] = D[:δt];

	return D
	
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
		val = .1*(s+t);
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = λ(newpt;case=:st);
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
		val = .5*s;
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = β(newpt;case=:st);
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
		val = exp(-s/50);
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = α(newpt;case=:st);
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
		
		# Defintion of α(s,t) given here	
		val = 2*s;
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = γ(newpt;case=:st);
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
		
		# Defintion of α(s,t) given here	
		val = 1.;
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = fˢ(newpt;case=:st);
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
		
		# Defintion of α(s,t) given here	
		val = 1.;
	
	elseif case == :χτ
		newpt = Fχτ(pt);
		val = fⁱ(newpt;case=:st);
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
		val[i] = λ(@view pts[:,i],prm;case=case);
	end

	return val
end
function β(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		val[i] = β(@view pts[:,i],prm;case=case);
	end

	return val
end
function α(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		val[i] = α(@view pts[:,i],prm;case=case);
	end

	return val
end
function γ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		val[i] = γ(@view pts[:,i],prm;case=case);
	end

	return val
end
function fˢ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		val[i] = fˢ(@view pts[:,i],prm;case=case);
	end

	return val
end
function fⁱ(pts::Matrix{Float64},prm::Dict{Symbol,Float64};case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		val[i] = fⁱ(@view pts[:,i],prm;case=case);
	end

	return val
end
