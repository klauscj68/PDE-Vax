# Collection of routines for solving the scalar conservation laws when the
# other two dependent variables are fixed.

# finterp
"""
Structure for interpolating and storing α,β,γ,λ on the rectangular mesh. 
These functions are defined in the inner constructor. The functions are
stored as a vector of nodal values ordered like the nodes of mesh.
"""
struct finterp 
	α::Vector{Float64}
	β::Vector{Float64}
	γ::Vector{Float64}
	λ::Vector{Float64}
	fˢ::Vector{Float64}
	fⁱ::Vector{Float64} # actually have be ρfⁱ?

	function finterp(mymesh::mesh)
		nnd = length(mymesh.nd);

		α = Vector{Float64}(undef,nnd);
		β = Vector{Float64}(undef,nnd);
		γ = Vector{Float64}(undef,nnd);
		λ = Vector{Float64}(undef,nnd);
		fˢ= Vector{Float64}(undef,nnd);
		fⁱ= Vector{Float64}(undef,nnd);

		for i=1:nnd
			s,t = mymesh.nd[i];
			α[i] = maximum([(365. -t)/365.,0.]);
			β[i] = .1*s;
			γ[i] = .05*t;
			λ[i] = .05*s*t;
			fˢ[i] = 1.;
			fⁱ[i] = 1.;
		end

		return new(α,β,γ,λ,fˢ,fⁱ)
	end
		
end

# eval
"""
An eval routine that evaluates the yˢ,yᵛ,yⁱ according to its conservation
law solution when the other two are fixed and boundary data is given.
sol:: Symbol matching :yˢ,:yᵛ,:yⁱ for the case being solved
∂s:: stores nodal values of s-axis ∂ data (matching a mesh)
∂t:: stores nodal values of t-axis ∂ data (matching a mesh)
mymesh:: mesh structure defining the geometry
ys:: dictionary storing the nodal values, matching mymesh.nd, for the two
     y's that are not solved for and if case is :yⁱ additional values at keys
     :yⁱ∂s and :yⁱ∂t for its boundary data which is solved elsewhere
frc:: all the α,β,γ,λ,fˢ,fⁱ terms are stored here
q:: optional arg, vector of nodes saying where to evaluate the solution. 
    If none are given it defaults to the mesh nodes
"""
function eval(sol::Symbol,mymesh::mesh,ys::Dict{Symbol,Vector{Float64}},
	      frc::finterp; q::Vector{node}=mymesh.nd)
	nq = length(q);
	nnd = length(mymesh.nd);
	gaussqd = quad1d();

	val = Vector{Float64}(undef,nq);

	if sol == :yˢ
		## Compute exponential argument's ∫'s
		# first integrand
		f₁du = frc.β.*ys[:yⁱ];

		# first interior integral at all nodes (s,t)
		∫f₁du = Vector{Float64}(undef,nnd);
		dir = [1.,0.];
		for i=1:nnd
			t = mymesh.nd[i].pt[2];
			pₒ = [0.,t];
			∫f₁du[i] = ∫fdτ(f₁du,mymesh,
					dir,mymesh.saxis,pₒ;
					gaussqd=gaussqd);
		end

		# second integrand
		f₂dv = frc.λ+∫f₁du;

		# second outside integral at all query points
		∫f₂dv = Vector{Float64}(undef,nq);
		dir = [1. /sqrt(2),1. /sqrt(2)];
		for i=1:nq
			pₒ = q[i].∂pt;
			# find last χaxis point underneath this δl value 
			# for dom discretization
			pos = 1;
			for j=1:mymesh.ntics
				if q[i].δl > mymesh.χaxis[j]
					pos += 1;
				else
					break
				end
			end
			χaxis = (q[i].δl == mymesh.χaxis[pos]) ? mymesh.χaxis : [mymesh.χaxis;q[i].δl];
			χaxis *= 1. /sqrt(2);
			
			# compute second integral
			∫f₂dv[i] = ∫fdτ(f₂dv,mymesh,
					dir,χaxis,pₒ;
					gaussqd=gaussqd);
		end

		# compute solution
		for i=1:nq
			val[i] = ( (q[i].∂pt[2] == 0.) ? eval(frc.fˢ,mymesh,q[:,i]) : 0. )*
			           exp(-∫f₂dv[i]);
		end

	elseif sol == :yᵛ
		## Compute exponential argument's ∫'s
		# first integrand
		f₁du = frc.β.*ys[:yⁱ];

		# first interior integral at all nodes (s,t)
		∫f₁du = Vector{Float64}(undef,nnd);
		dir = [1.,0.];
		for i=1:nnd
			t = mymesh.nd[i].pt[2];
			pₒ = [0.,t];
			∫f₁du[i] = ∫fdτ(f₁du,mymesh,
					dir,mymesh.saxis,pₒ;
					gaussqd=gaussqd);
		end

		# second integrand
		f₂dv = frc.α + (1-frc.α)*∫f₁du;

		# second outside integral at all query points
		∫f₂dv = Vector{Float64}(undef,nq);
		dir = [1. /sqrt(2),1. /sqrt(2)];
		for i=1:nq
			pₒ = q[i].∂pt;
			# find last χaxis point underneath this δl value 
			# for dom discretization
			pos = 1;
			for j=1:mymesh.ntics
				if q[i].δl > mymesh.χaxis[j]
					pos += 1;
				else
					break
				end
			end
			χaxis = (q[i].δl == mymesh.χaxis[pos]) ? mymesh.χaxis : [mymesh.χaxis;q[i].δl];
			χaxis *= 1. /sqrt(2);
			
			# compute second integral
			∫f₂dv[i] = ∫fdτ(f₂dv,mymesh,
					dir,χaxis,pₒ;
					gaussqd=gaussqd);
		end

		# compute solution
		dir = [1.,0];
		for i=1:nq
			fdτ = frc.λ.*ys[:yˢ];
			pₒ = [0.,q[i].pt[2]];
			y∂t = ∫fdτ(fdτ,mymesh,dir,mymesh.saxis,pₒ;
				   gaussqd=gaussqd);

			val[i] = ( (q[i].∂pt[2] == 0.) ? 0. : y∂t )*
			           exp(-∫f₂dv[i]);
		end

	elseif sol == :yⁱ
		dir = [1. /sqrt(2),1. /sqrt(2)];
		for i=1:nq
			pₒ = q[i].∂pt;
			# find last χaxis point underneath this δl value
			# # for dom discretization
			pos = 1;
			for j=1:mymesh.ntics
				if q[i].δl > mymesh.χaxis[j]
					pos += 1;
				else
					break
				end
				χaxis = (q[i].δl == mymesh.χaxis[pos]) ? mymesh.χaxis : [mymesh.χaxis;q[i].δl];
				χaxis *= 1. /sqrt(2);
			end
			# Compute the integral in the exponential argument
			∫fdv = ∫fdτ(frc.γ,mymesh,dir,χaxis,pₒ;
			            gaussqd=gaussqd);
			y∂p = (q[i].∂pt[1] == 0.) ? myinterp(mymesh.taxis,ys[:yⁱ∂t],q[i].∂pt[2]) : myinterp(mymesh.saxis,ys[:yⁱ∂s],q[i].∂pt[1]);
			
			val[i] = y∂p*exp(∫fdv);
		end

	else
		error("Not a valid solution case")
	end
	
	return val

end

# yspl
"""
Structure for storing solution splines via their boundary data
sol:: is either :yˢ,:yᵛ,:yⁱ depending on the case
∂s:: stores nodal values of s-axis boundary values (matching a mesh)
∂t:: stores nodal values of t-axis boudary values (matching a mesh)
c:: stores nodal values across the full rectangular mesh

ys:: Dict storing the nodal values of the other two ysol's which are 
     fixed for the purpose of solving for the one specified by sol
"""
struct yspl
	sol::Symbol
	∂s::Vector{Float64}
	∂t::Vector{Float64}
	c::Vector{Float64}
	
	#function yspl(sol::Symbol,mymesh::mesh,
	#	      ys::Dict{Symbol,Vector{Float64}}

end
