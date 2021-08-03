# Collection of routines for optimizing the boundary conditions to satisfy
# recursive system

#%% Errors
# ∂err
"""
Compute the ℓ₂² error between prescribed boudary conditions and values of 
implicit conditions they are supposed to satisfy.
"""
function ∂err(y0::Dict{Symbol,Vector{Float64}},YSOL::Dict{Symbol,Vector{Vector{Float64}}},
		dom::Domain,prm::Dict{Symbol,Float64})
	# Number of nonzero χaxis points
	nnd = dom.nelm+1;
	nτs = length(YSOL[:hτs][1]);
	∂χ = @view dom.χaxis[1:nnd];

	#  Initialize Gbpts and gaussqd
	gaussqd = quad1d(prm[:nqd]);
	Gbpts = Vector{Matrix{Float64}}(undef,prm[:nelm]);
	for i=1:prm[:nelm]
		Gbpts[i] = Matrix{Float64}(undef,2,prm[:nqd]);
	end

	# Compute remaining ∫'s all at (0,t) in (s,t)
	λyˢds = Vector{Vector{Float64}}(undef,nτs);
	Imαyᵛds = Vector{Vector{Float64}}(undef,nτs);

	@inbounds for i=1:nτs
		λyˢds[i] = Vector{Float64}(undef,dom.nnd);
		Imayᵛds[i] = Vector{Float64}(undef,dom.nnd);
	end

	λyˢ!(λyˢds,YSOL[:yˢ],YSOL[:hτs][1],dom)
	Imayᵛ!(Imayᵛds,YSOL[:yᵛ],YSOL[:hτs][1],dom);

	
	# compute error
	ℓ₂² = 0.;
	χτ = [0.,0.];
	@inbounds for i=1:nnd
		χτ[1] = ∂χ[i];
		
		# yᵛ error
		ℓ₂² += ( pullb∫fds!(λyˢds,YSOL[:hτs][1],χτ,dom,prm[:nelm],Gbpts;
				    gaussqd=gaussqd) -
			      y0[:yᵛ][i] )^2;

		# yⁱ error 
		∫yˢ = pullb∫fds!(YSOL[:yˢ],YSOL[:hτs][1],χτ,dom,prm[:nelm],Gbpts;gaussqd=gaussqd);
		∫βy = YSOL[:∫βyⁱds][1][i]; #τ=0 point
		∫Imαyᵛ = pullb∫fds!(Imayᵛds,YSOL[:hτs][1],χτ,dom,prm[:nelm],Gbpts;gaussqd=gaussqd);

		ℓ₂² += ( (∫yˢ+∫Imαyᵛ)*∫βy - y0[:yⁱ][i] )^2;
	end
	
	return ℓ₂²
end
