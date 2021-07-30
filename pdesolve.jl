# Routines used to integrate the vaccination pde system without coupled ∂-data

#%% ODE solver
#flowfield!
"""
Evaluate the flowfield derivative at a point all in (χ,τ) coordinate plane
Mutates the ∂τy and Gbpts values. Note that y is a vertically stacked [yˢ,yᵛ,yⁱ] 
vector of the nodals.

∂τy:: in-place derivative value
y:: Value of solution at point χτ
H:: A dictionary storing the past history values in (χ,τ) coordinates of y,
    β at the nodal values for sequential timesteps hτs
    Keys are :y,:β. Mutates the value of H[:y][nτstep] to y
hτs:: The τ-time values at which history is stored. hτs[1:nτstep-1] is 
      histories solved for already. Mutates the value of hτs[nτstep] to now
nτstep:: The index of integration downsample step we are NOW solving for
         (only changes when accumulated enough integration steps to pass
	  a downsample)
dom:: The domain over which integration
Gbpts:: in-place storage for guassian quadrature points
βyⁱ:: in-place storage size [>=nτstep][nnd] for computing product terms
      β*yⁱ 
      Note: be better if you could allocate whole downsample storage before 
            integrating whole length and use slices
"""
function flowfield!(∂τy::VecVw,y::VecVw,τ::Float64,
		    H::Dict{Symbol,VecVecVw},hτs::VecVw,nτstep::Int64, # The H dictionary has wrong type
		    prm::Dict{Symbol,Float64},
		    dom::Domain,
		    Gbpts::Vector{Matrix{Float64}},
		    βyⁱ::Vector{Vector{Float64}},
		    gaussqd::Dict{Symbol,Vector{Float64}}=quad1d())
	# Initalize
	#  Update solution history to implicitly include the value being 
	#  solved for
	H[:y][nτstep] = y; hτs[nτstep] = τ;
	nelm = prm[:nelm];

	yˢ = @view y[1:dom.nnd];
	yᵛ = @view y[dom.nnd+1:2*dom.nnd];
	yⁱ = @view y[2*dom.nnd+1:3*dom.nnd];

	∂τyˢ = @view ∂τy[1:dom.nnd];
	∂τyᵛ = @view ∂τy[dom.nnd+1:2*dom.nnd];
	∂τyⁱ = @view ∂τy[2*dom.nnd+1:3*dom.nnd];
	
	# Compute βyⁱ 
	@inbounds for i=1:nτstep
		βyⁱ[i] = (@view H[:β][i]).*
		           (@view H[:y][i][2*dom.nnd+1:3*dom.nnd]);
	end

	# Compute pullback line integral ∫ᴸ₀βyⁱds 
	# Temp store in-place to ∂τyⁱ
	∫ᴸ₀βyⁱds = ∂τyⁱ;
	@inbounds for i=1:dom.nnd
		∫ᴸ₀βyⁱds[i] = pullb∫fds!(@view βyⁱ[1:nτstep],
					 @view hτs[1:nτstep],[dom.χaxis[i],τ],dom,
					 nelm,Gbpts;gaussqd=gaussqd);
	end
	
	# Compute derivatives
	pt = [0.,0.];
	@inbounds for i=1:dom.nnd
		pt[:] = [dom.χaxis[i],τ];
		αs = α(pt;case=:χτ);

		∂τyˢ[i] = -(λ(pt;case=:χτ) + ∫ᴸ₀βyⁱds[i])*yˢ[i];	
		∂τyᵛ[i] = -(αs + (1-αs)*∫ᴸ₀βyⁱds[i])*yᵛ[i];
		∂τyⁱ[i] = -γ(pt;case=:χτ);
	end

	∂τy *= 1/sqrt(2);
end

# runga!
""" 
6th order explicit runga-kutta solver to use for numerical integration
"""
function runga!(K0::VecWv,K1::VecWv,K2::VecWv,K3::VecWv,
		∂τy::VecVw,y::VecVw,τ::Float64,δτ
		    H::Dict{Symbol,VecVecVw},hτs::VecVw,nτstep::Int64, # The H dictionary has wrong type
		    prm::Dict{Symbol,Float64},
		    dom::Domain,
		    Gbpts::Vector{Matrix{Float64}},
		    βyⁱ::Vector{Vector{Float64}},
		    gaussqd::Dict{Symbol,Vector{Float64}})

        K0[:] = δτ*flowfield(sheet,mydep,τ,y,frc_M);
        K1[:] = δτ*flowfield(sheet,mydep,τ+0.5*δτ,y+0.5*K0,frc_M);
        K2[:] = δτ*flowfield(sheet,mydep,τ+0.5*δτ,y+0.5*K1,frc_M);
        K3[:] = δτ*flowfield(sheet,mydep,τ+δτ,y+K2,frc_M)
end

# odesolve
"""
An ODE integrator for the vaccination PDE system before coupling ∂-data
"""
function odesolve(dom::Domain,∂y::Dict{Symbol,VecVw},prm::Dict{Symbol,Float64}) # ∂y has wrong type
	# Initialize
	#  Initial data
	y0 = [∂y[:yˢ];
	      ∂y[:yᵛ];
	      ∂y[:yⁱ]];
	
	#  τ-downsample axis
	τdwn = convert(Vector,0.: prm[:δτ] : prm[:τfin]);

	#  Preallocate for history and store y0,β
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	H = Dict{Symbol,VecVecVw}();
	#    y0
	H[:y] = Vector{Vector{Float64}}(undef,length(τdwn)); # has wrong type MAYBE THIS TYPE ENOUGH
	@inbounds for i=1:length(τdwn)
		H[:y][i] = Vector{Float64}(undef,3*dom.nnd);
	end 
	H[:y][1] = y0;
	
	#    β
	H[:β] = Vector{Vector{Float64}}(undef,length(τdwn)); 
	@inbounds for i=1:length(τdwn)
		H[:β][i] = Vector{Float64}(undef,dom.nnd);
	end
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	
	gen = [1,1];
	@inbounds for k=1:length(τdwn)*dom.nnd
		# Take a view of things adjacent in memory
		i = gen[1]; j = gen[2];
		if j == 1
			ram = @view H[:β][i]
		end
		
		# Compute β at point and store
		pt = [dom.χaxis[j],τdwn[i]];
		ram[j] = β(pt;case=:χτ);

		# cycle generator: i τdwnsmp and j is nodal index
		if j != dom.nnd
			gen[2] += 1;
		else
			gen[1] += 1;
			gen[2] = 1;
		end
	end
	
	#  Preallocate for gaussqd and Gbpts
	nelm = Int64(prm[:nelm]); nqd = Int64(prm[:nqd]);
	gaussqd = quad1d(nqd);
	Gbpts = Vector{Matrix{Float64}}(undef,nelm);
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	@inbounds for i=1:nelm
		Gbpts[i] = Matrix{Float64}(undef,2,nqd);
	end
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	# preallocate for βyⁱ
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	βyⁱ = Vector{Vector{Float64}}(undef,length(τdwn));
	for i=1:length(τdwn)
		βyⁱ[i] = Vector{Float64}(undef,dom.nnd);
	end
	#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
