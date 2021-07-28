# Routines used to integrate the vaccination pde system without coupled ∂-data

#%% ODE solver
#flowfield!
"""
Evaluate the flowfield derivative at a point all in (χ,τ) coordinate plane
Mutates the ∂y and Gbpts values. Note that y is a vertically stacked [yˢ,yᵛ,yⁱ] 
vector of the nodals.

∂y:: in-place derivative value
y:: Value of solution at point χτ
H:: A dictionary storing the past history values in (χ,τ) coordinates of y,
    β at the nodal values for sequential timesteps hτs
    Keys are :y,:β
hτs:: The τ-time values at which history is stored
τ:: The integration time at which to evaluate the solution
dom:: The domain over which integration
Gbpts:: in-place storage for guassian quadrature points
βyⁱ:: in-place storage size [>=length(hτs)][nnd] for computing product terms
      β*yⁱ 
      Note: be better if you could allocate whole downsample storage before 
            integrating whole length and use slices
"""
function flowfield!(∂y::VecVw,y::VecVw,H::Dict{Symbol,VecVecVw},hτs::VecVw,
		    prm::Dict{Symbol,Float64},τ::Float64,
		    dom::Domain,nelm::Int64,
		    Gbpts::Vector{Matrix{Float64}},
		    βyⁱ::Vector{Vector{Float64}};
		    gaussqd::Dict{Symbol,Vector{Float64}}=quad1d())
	# Initalize
	nelm = length(Gbpts);

	yˢ = @view y[1:dom.nnd];
	yⁱ = @view y[dom.nnd+1:2*dom.nnd];
	yᵛ = @view y[2*dom.nnd+1:3*dom.nnd];

	∂yˢ = @view ∂y[1:dom.nnd];
	∂yⁱ = @view ∂y[dom.nnd+1:2*dom.nnd];
	∂yᵛ = @view ∂y[2*dom.nnd+1:3*dom.nnd];
	
	# Compute βyⁱ
	@inbounds for i=1:length(hτs)
		βyⁱ[i] = (@view H[:β][i]).*
		           (@view H[:y][i][dom.nnd+1:2*dom.nnd]);
	end

	# Compute pullback line integral ∫ᴸ₀βyⁱds 
	# Temp store in-place to ∂yⁱ
	∫ᴸ₀βyⁱds = ∂yⁱ;
	@inbounds for i=1:dom.nnd
		∫ᴸ₀βyⁱds[i] = pullb∫fds!(@view βyⁱ[1:length(hτs)],
					 hτs,[dom.χaxis[i],τ],dom,
					 nelm,Gbpts;gaussqd=gaussqd);
	end
	
	# Compute derivatives
	pt = [0.,0.];
	@inbounds for i=1:dom.nnd
		pt[:] = [dom.χaxis[i],τ];
		αs = α(pt;case=:χτ);

		∂yˢ[i] = -(λ(pt;case=:χτ) + ∫ᴸ₀βyⁱds[i])*yˢ[i];	
		∂yᵛ[i] = -(αs + (1-αs)*∫ᴸ₀βyⁱds[i])*yᵛ[i];
		∂yⁱ[i] = -γ(pt;case=:χτ);
	end

	∂y *= 1/sqrt(2);
end
