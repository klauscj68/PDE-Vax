## File storing parameters and model equation terms
# Should also have loaded auxilliary.jl

#%% Model scalar parameters
# data
"""
Dictionary storing model scalar parameters
"""
function data()
	prm = DSymFl();
	
	# Domain
	#  Largest age by yˢ,yᵛ,yⁱ
	prm[:Ls] = 100. *365; # up to 100 years
	prm[:Lv] = 11.; # up to 8 months
	prm[:Li] = 45.; # up to 21 days

	#  Largest time
	prm[:T] = 31.;	

	# Epidemic
	#  Initial fraction vaccinated
	prm[:ρ] = .3;

	# Hazard rates
	#  β: Weibull distribution
	prm[:βa] = .5;
	prm[:βb] = 10.;
	prm[:βη] = 1.; # mutated by ∂YSOL! to ensure continuity at (0,0)

	#  γ: Weibull distributrion
	prm[:γa] = 20.;
	prm[:γb] = 1.5;

	#  Initial populations
	prm[:fˢη] = 1.; # mutated by ∂YSOL! to ensure prob distribution
	prm[:fⁱη] = 1.; # mutated by ∂YSOL! to ensure prob distribution

	# Numerical discretization
	#  Number nodes within each [t=t₀] set
	prm[:nnd] = 1225.;

	#  downsamples used to save the solution along integration in space-time
	prm[:δt] = .1;
	prm[:δs] = .01; # relative length of axis

	#   tolerances for ode integration and maximum permitted 
	#   integration step
	prm[:atol] = 1e-9;
	prm[:rtol] = 1e-3;
	prm[:rlow] = 1e-6; # least yval for rel error exactly enforced
	prm[:δtmax] = 3.;

	return prm
	
end
# data!
"""
Mutate the prm dictionary so that all parameters depending on others are 
properly set to those values
"""
function data!(prm::DSymFl)
	#-----
	# Compute normalization constants
	nnd = Int64(prm[:nnd]);

	#  Define initial geometry
	tlvl_yˢ = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:Ls],nnd)) );
	tlvl_yᵛ = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:Lv],nnd)) );
	tlvl_yⁱ = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:Li],nnd)) );

	vβ = Vector{Float64}(undef,nnd);
	vfˢ = Vector{Float64}(undef,nnd);
	vfⁱ = Vector{Float64}(undef,nnd);
	
	#  Adjust fˢ,fⁱ to be probability distributions
	@inbounds for i=1:nnd
		nd = @view tlvl_yˢ.nds[:,i];
		vfˢ[i] = fˢ(nd,prm;case=:χτ);

		nd = @view tlvl_yⁱ.nds[:,i];
		vfⁱ[i] = fⁱ(nd,prm;case=:χτ);
	end
	
	#   Record the scaling
	∫fˢds = ∫line(Yℓvℓ(tlvl_yˢ,vfˢ));
	∫fⁱds = ∫line(Yℓvℓ(tlvl_yⁱ,vfⁱ));

	prm[:fˢη] *= 1/∫fˢds;
	prm[:fⁱη] *= 1/∫fⁱds; vfⁱ *= prm[:fⁱη];
	
	#  Adjust the β to be compatible with cont BC's at (0,0)
	@inbounds for i=1:nnd
		nd = @view tlvl_yⁱ.nds[:,i];
		vβ[i] = β(nd,prm;case=:χτ);
	end

	∫βfⁱds = ∫line(Yℓvℓ(tlvl_yⁱ,vβ.*vfⁱ));
	prm[:βη] *= fⁱ([0.,0.],prm;case=:χτ)/∫βfⁱds;
	prm[:βη] = (!isnan(prm[:βη])) ? (prm[:βη]) : 0.;

end

#%% Model function parameters
# λ
"""
Evaluate λ 
Defaults to λ(s,t) but optional case argument can be used to instead
compute λ(s(χ,τ),t(χ,τ))
case:: can be either :st or :χτ
"""
function λ(pt::VecVw,prm::DSymFl;case::Symbol=:st)
	
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of λ(s,t) given here
		#  Note: need λ(s,0)≡0 for BC's
		val = 2 - abs(t/31 - 2);
		val = val >= 0 ? val : 0.;
	
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
function β(pt::VecVw,prm::DSymFl;case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of β(s,t) given here
		ram = s/prm[:βa];
		val = (ram >= 1e-6) ? prm[:βb]/prm[:βa]*real( (ram+0im)^(prm[:βb]-1.) ) : (
						prm[:βb]/prm[:βa]*real( (1e-6+0im)^(prm[:βb]-1.) ) );
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
function α(pt::VecVw,prm::DSymFl;case::Symbol=:st)
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
function γ(pt::VecVw,prm::DSymFl;case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of γ(s,t) given here
		ram = s/prm[:γa];
		val = (ram >= 1e-6) ? prm[:γb]/prm[:γa]*real( (ram+0im)^(prm[:γb]-1.) ) : (
			                 prm[:γb]/prm[:γa]*real( (1e-6+0im)^(prm[:γb]-1.) ) );
	
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
function fˢ(pt::VecVw,prm::DSymFl;case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of fˢ given here
		#  Match to Ohio age distribution (yrs)
		#  ∂YSOL! will compute correct normalization to take
		#  density from years -> days
		age = s/365;
		nd = [0.,5.,15.,25.,35.,45.,55.,65.,100.];
		ρ = [0.,.0125,.0125,.0125,.0125,.0125,.0125,.0125,0.];
		val = myinterp(nd,ρ,age);
		val *= prm[:fˢη];

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
function fⁱ(pt::VecVw,prm::DSymFl;case::Symbol=:st)
	if case == :st
		s = pt[1]; t = pt[2];
		
		# Defintion of fⁱ given here	
		val = 1+5-abs(s-5);
		val = val >= 0 ? val : 0.;
		val *= prm[:fⁱη];

	elseif case == :χτ
		newpt = Fχτ(pt);
		val = fⁱ(newpt,prm;case=:st);
	else
		error("not valid fⁱ-eval case");
	end

	return val
end

#%% Auxilliary methods for sampling function parameters across batch of sample points
function λ(pts::Matrix{Float64},prm::DSymFl;case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = λ(pt,prm;case=case);
	end

	return val
end
function β(pts::Matrix{Float64},prm::DSymFl;case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = β(pt,prm;case=case);
	end

	return val
end
function α(pts::Matrix{Float64},prm::DSymFl;case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = α(pt,prm;case=case);
	end

	return val
end
function γ(pts::Matrix{Float64},prm::DSymFl;case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = γ(pt,prm;case=case);
	end

	return val
end
function fˢ(pts::Matrix{Float64},prm::DSymFl;case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = fˢ(pt,prm;case=case);
	end

	return val
end
function fⁱ(pts::Matrix{Float64},prm::DSymFl;case::Symbol=:st)
	@assert size(pts)[1] == 2 "pts must have dimension 2"
	npts = size(pts)[2];

	val = Vector{Float64}(undef,npts);
	@inbounds for i=1:npts
		pt = @view pts[:,i];
		val[i] = fⁱ(pt,prm;case=case);
	end

	return val
end
