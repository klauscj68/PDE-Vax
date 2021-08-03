# Routines used to integrate the vaccination pde system without coupled ∂-data
using DifferentialEquations

#%% ODE solver
# flow∂yⁱ!
"""
Compute the flowfield equation for yⁱ
"""
function flow∂yⁱ!(∂yⁱ::Vector{Float64},yⁱ::Vector{Float64},p::Dict{Symbol,Vector{Vector{Float64}}},τ::Float64)
	χaxis = p[:χaxis];	
	@inbounds for i=1:length(∂yⁱ)
		pt = [χaxis[i],τ];
		∂yⁱ[i] = -gamma(pt;case=:χτ)*yⁱ[i];
	end
end

# flow∂yˢᵛ!
"""
Compute the flowfield equation for yˢ,yᵛ
"""
function flow∂yˢᵛ!(∂yˢᵛ::Vector{Float64},yˢᵛ::Vector{Float64},p::Dict{Symbol,Vector{Vector{Float64}}},τ::Float64)	

	# Initialize
	χaxis = p[:χaxis][1];
	nnd = length(χaxis);
	hτs = p[:hτs][1];
	∫βyⁱds = p[:∫βyⁱds];


	yˢ = @view yˢᵛ[1:nnd];
	yᵛ = @view yˢᵛ[nnd+1:2*nnd];

	∂yˢ = @view ∂yˢᵛ[1:nnd];
	∂yᵛ = @view ∂yˢᵛ[nnd+1:2*nnd];

	# Interpolate ∫βyⁱds @ time τ by in-place assignment to ∂yˢᵛ
	if τ <= hτs[1]
		∂yˢ[:] = ∫βyⁱds[1];
	elseif τ >= hτs[end]
		∂yˢ[:] = ∫βyⁱds[end];
	else
		pos = myfindfirst(hτs,τ);
		η = (τ - hτs[pos-1])/(hτs[pos]-hτs[pos-1]);
		∂yˢ[:] = (1-η)*∫βyⁱds[pos-1] .+ η*∫βyⁱds[pos]; 
	end

	∂yᵛ[:] = ∂yˢ;

	# Compute system
	χτ = [0.,τ];
	@inbounds for i=1:nnd
		χτ[1] = χaxis[i];
		# yˢ contribution
		∂yˢ[i] += λ(χτ;case=:χτ); ∂yˢ[i] *= -1; ∂yˢ[i] *= yˢ[i];

		# yᵛ contribution
		α = α(χτ;case=:χτ);
		∂yᵛ[i] *= (1-α); ∂yᵛ[i] += α; ∂yᵛ[i] *= -1; ∂yᵛ[i] *= yᵛ[i];

	end
	
end

# modelrun
"""
Run the model before coupling ∂-data. y0 is the ∂-data
"""
function modelrun(y0::Dict{Symbol,Vector{Float64}},prm::Dict{Symbol,Float64},dom::Domain)
	τspan = [0.,prm[:τfin]];

	# Solve the yⁱ system
	p0 = Dict{Symbol,Vector{Vector{Float64}}}();
	p0[:χaxis] = [dom.χaxis];

	prob_yⁱ = ODEProblem(flow∂yⁱ!,y0[:yⁱ],τspan,p0);
	sol = solve(prob_yⁱ,Tsit5(),reltol=prm[:rtol],abstol=prm[:atol],saveat=prm[:δτ]);
	
	yⁱ = sol.u;
	hτs = [sol.t]; nτs = length(hτs[1]);

	# Store solution
	YSOL = Dict{Symbol,Vector{Vector{Float64}}}(:yⁱ=>yⁱ,:hτs=>hτs);

	# Compute the βyⁱ term
	βyⁱ = Vector{Vector{Float64}}(undef,nτs);
	for i=1:dom.nnd
		βyⁱ[i] = Vector{Float64}(undef,dom.nnd);
	end
	βyⁱ!(βyⁱ,yⁱ,hτs[1],dom);

	# Compute ∫βyⁱds
	#  Initialize ∫βyⁱds
	∫βyⁱ = deepcopy(βyⁱ);

	#  Initialize Gbpts and gaussqd
	gaussqd = quad1d(prm[:nqd]);
	Gbpts = Vector{Matrix{Float64}}(undef,prm[:nelm]);
	for i=1:prm[:nelm]
		Gbpts[i] = Matrix{Float64}(undef,2,prm[:nqd]);
	end

	χτ = [0.,0.];
	gen = [1,1];
	for k=1:nτs*dom.nnd
		# i is τ-step and j is the dom.nnd
		i = gen[1]; j = gen[2];
		if j == 1
			∫βnow = ∫βyⁱ[i];
			χτ[2] = hτs[i];
		end

		χτ[1] = dom.χaxis[j];

		# Compute pullback value
		∫βnow[j] = pullb∫fds!(βyⁱ,hτs[1],χτ,dom,prm[:nelm],Gbpts;gaussqd);

		# cycle generator
		if j != dom.nnd
			gen[2] += 1;
		else
			gen[1] += 1;
			gen[2] = 1;
		end
	end

	# Compute solution to yˢ,yᵛ
	p∞ = Dict{Symbol,Vector{Vector{Float64}}}(:hτs=>hτs,:∫βyⁱds=>∫βyⁱ,:χaxis=>[dom.χaxis]);
	prob_yˢᵛ = ODEProblem(flow∂yˢᵛ!,[y0[:yˢ];y0[:yᵛ]],τspan,p∞);
	sol∞ = solve(prob_yˢᵛ,Tsit5(),reltol=prm[:rtol],abstol=prm[:atol],saveat=prm[:δτ]);

	yˢᵛ = sol∞.u;
	yˢ = Vector{Vector{Float64}}(undef,nτs);
	yᵛ = Vector{Vector{Float64}}(undef,nτs);

	for i=1:nτs
		yˢ[i] = yˢᵛ[i][1:dom.nnd];
		yᵛ[i] = yˢᵛ[i][dom.nnd+1:2*dom.nnd];
	end
	
	YSOL[:yˢ] = yˢ;
	YSOL[:yᵛ] = yᵛ;
	YSOL[:∫βyⁱds] = ∫βyⁱ;
	
	return YSOL
end
