# Collection of routines for solving the Fredholm ∂-condition of yⁱ.
#  Requires mesh.jl, splines.jl, CLscalar.jl

# Γ
"""
Ancillary function used to compute the Fredholm kernel
Γ(t) = ∫∞₀[ys(s,t) + (1-α(s))*yᵛ(s,t)]ds

returns Γ evaluated at the mesh nodes
"""
function Γ(mymesh::mesh,
	   ys::Dict{Symbol,Vector{Float64}},frc::finterp;
	   gaussqd::Dict{Symbol,Vector{Float64}}=quad1d())

	fds = ys[:yˢ]+(1-frc.α).*ys[:yᵛ];
	
	∫fds = Vector{Float64}(undef,length(mymesh.taxis));
	dir = [1.,0.];
	for i=1:length(mymesh.taxis)
		∫fds[i] = ∫fdτ(fds,mymesh,dir,mymesh.saxis,
			       [0.,mymesh.taxis[i]];
			       gaussqd);
	end

	# Reshape into nodal value vector
	val = repeat(reshape(∫fds,(1,length(mymesh.taxis))),
			     outer=(length(mymesh.saxis),1));
	return reshape(val,length(mymesh.nd));
end

#Σ
"""
Ancillary function used to return the Fredholm forcing term
Σ(t) = (∫^∞_t β(u)*ρ*fⁱ(u)exp(-∫ᵗ₀γ(u+v)dv)du)
"""
function Σ(mymesh::mesh,
	   ys::Dict{Symbol,Vector{Float64}},frc::finterp;
	   gaussqd::Dict{Symbol,Vector{Float64}}=quad1d())
	
	nnd = length(mymesh.nd);

	f2dv = frc.γ;
	∫f2dv = Vector{Float64}(undef,nnd)
	dir = [0.,1.];
	for i=1:nnd
		pₒ = mymesh.nd[i];
		dom = ∫dom(mymesh.taxis,[0.,pₒ[2]]);
		∫f2dv = ∫fdτ(f2dv,mymesh,dir,dom,pₒ;
			     gaussqd);
	end

	f1du = frc.β.*frc.fⁱ.*exp.(-∫f2dv);

	∫f1du = Vector{Float64}(undef,length(mymesh.saxis));
	dir = [1.,0.];
	for i=1:length(mymesh.taxis)
		pₒ = mymesh.nd[i];
		dom = ∫dom(mymesh.saxis,[pₒ[2],mymesh.urg[2]]);
		∫f1du[i] = ∫fdτ(f1du,mymesh,dir,dom,pₒ,
				  gaussqd);
	end

	val = repeat(reshape(∫f1du,(1,length(mymesh.taxis))),
			     outer=(length(mymesh.saxis),1));
	return reshape(val,length(mymesh.nd));

end


# FrdK
"""
Compute the value of the Fredholm kernel at all mesh nodes whose value is
K(s,t) = Γ(t)
"""

