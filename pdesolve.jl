# Routines used to integrate the vaccination pde system
using CSV,DataFrames

#%% Ancillary routines for computing system integrands
# βy!
"""
Compute βy at the nodes along Tℓvℓ. Routine mutates the optional
argument βy.
"""
function βy!(ylvl::Yℓvℓ,prm::DSymFl;
	     βy::Yℓvℓ=Yℓvℓ(undef))

	if isnan(βy.ys[1])
		βy = deepcopy(ylvl);
		flagrt = true;
	else
		flagrt = false;
	end
	
	for i=1:ylvl.tlvl.nnd
		nd = @view ylvl.tlvl.nds[:,i];
		βy.ys[i] = β(nd,prm;case=:χτ)*
		             ylvl.ys[i];
	end

	if flagrt
		return βy
	end
end

# λy!
"""
Compute λy at the nodes along Tℓvℓ. Routine mutates the optional 
argument λy
"""
function λy!(ylvl::Yℓvℓ,prm::DSymFl;
	     λy::Yℓvℓ=Yℓvℓ(undef))

	if isnan(λy.ys[1])
		λy = deepcopy(ylvl);
		flagrt = true;
	else
		flagrt = false;
	end

	for i=1:ylvl.tlvl.nnd
		nd = @view ylvl.tlvl.nds[:,i];
		λy.ys[i] = λ(nd,prm;case=:χτ)*
		             ylvl.ys[i];
	end

	if flagrt
		return λy
	end
end
		

# Imαy!
"""
Compute (1-α)y at the nodes along Tℓvℓ. Routine mutates the optional
argument Imαy
"""
function Imαy!(ylvl::Yℓvℓ,prm::DSymFl;
	       Imαy::Yℓvℓ=Yℓvℓ(undef))
	
	if isnan(Imαy.ys[1])
		Imαy = deepcopy(ylvl);
		flagrt = true;
	else
		flagrt = false;
	end

	for i=1:ylvl.tlvl.nnd
		nd = @view ylvl.tlvl.nds[:,i];
		Imαy.ys[i] = (1-α(nd,prm;case=:χτ))*
			      ylvl.ys[i];
	end

	if flagrt
		return Imαy
	end
end

#%% System integrator
# euler!
"""
Compute an explicit Euler step of the PDE vaccination system, mutates the
optional dictionary input for writing solutions and λ,α,γ and 
βy,λy,(1-α)y values

Note: Keep the EYSOL ∫yds fields NaN and then euler! will internally 
      mutate those integrals to their correct values when in later iteration
      it is passed as YSOL
"""
function euler!(δt::Float64,YSOL::Dict{Symbol,Yℓvℓ},prm::DSymFl;
		EYSOL::Dict{Symbol,Yℓvℓ}=Dict{Symbol,Yℓvℓ}(),
		updTℓvℓ::Vector{Symbol}=[:yˢ,:yᵛ,:yⁱ,:λ,:α,:γ,:βyⁱ,:λyˢ,:Imαyᵛ])
	
	nnd = YSOL[:yˢ].tlvl.nnd;
	t₀ = YSOL[:yˢ].tlvl.t₀[1];	
	if isempty(EYSOL)
		EYSOL = deepcopy(YSOL);
		flagrt = true;
	else
		flagrt = false;
	end
	
	# Compute needed integrals for forward Euler step
	#  Note the integral value is constant within [t=t₀]
	if isnan(YSOL[:βyⁱ].∫yds[1])
		∫βyⁱds = ∫line(YSOL[:βyⁱ]);
		YSOL[:βyⁱ].∫yds[1] = ∫βyⁱds;
	else
		∫βyⁱds = YSOL[:βyⁱ].∫yds[1];
	end

	if isnan(YSOL[:λyˢ].∫yds[1])
		∫λyˢds = ∫line(YSOL[:λyˢ]);
		YSOL[:λyˢ].∫yds[1] = ∫λyˢds;
	else
		∫λyˢds = YSOL[:λyˢ].∫yds[1];
	end

	if isnan(YSOL[:yˢ].∫yds[1])
		∫yˢds = ∫line(YSOL[:yˢ]);
		YSOL[:yˢ].∫yds[1] = ∫yˢds;
	else
		∫yˢds = YSOL[:yˢ].∫yds[1];
	end

	if isnan(YSOL[:Imαyᵛ].∫yds[1])
		∫Imαyᵛds = ∫line(YSOL[:Imαyᵛ]);
		YSOL[:Imαyᵛ].∫yds[1] = ∫Imαyᵛds;
	else
		∫Imαyᵛds = YSOL[:Imαyᵛ].∫yds[1];
	end

	# Advance to new tlvl's for this Euler step
	Tℓvℓ!(t₀+δt,YSOL[:yˢ].tlvl,EYSOL[:yˢ].tlvl);
	Tℓvℓ!(t₀+δt,YSOL[:yᵛ].tlvl,EYSOL[:yᵛ].tlvl);
	Tℓvℓ!(t₀+δt,YSOL[:yⁱ].tlvl,EYSOL[:yⁱ].tlvl);
	for key in updTℓvℓ
		if (key == :yˢ)||(key == :yᵛ)||(key == :yⁱ)
			continue
		elseif (key == :λ)||(key == :λyˢ)
			ytlvl = EYSOL[:yˢ].tlvl;
		elseif (key == :α)||(key == :Imαyᵛ)
			ytlvl = EYSOL[:yᵛ].tlvl;
		else
			ytlvl = EYSOL[:yⁱ].tlvl;
		end
		EYSOL[key].tlvl.t₀[:] = ytlvl.t₀;
		EYSOL[key].tlvl.nds[:,:] = ytlvl.nds;
		EYSOL[key].tlvl.χrg[:] = ytlvl.χrg;
		EYSOL[key].tlvl.τrg[:] = ytlvl.τrg;
	end
	
	# Update the ∂-nodal solution value at (χ,τ) = (-t₀,0)
	EYSOL[:yˢ].ys[1] = 0.;
	EYSOL[:yᵛ].ys[1] = ∫λyˢds;
	EYSOL[:yⁱ].ys[1] = (∫yˢds + ∫Imαyᵛds)*∫βyⁱds;
	
	# Update the non-∂ nodal solution values for this Euler step 
	gen = [1,1]; ykey = :yˢ; χsY = @view YSOL[ykey].tlvl.nds[1,:];
	@inbounds for i=1:3*(nnd-1)
		# Cycle generator gen[1] is for yˢ,yᵛ,yⁱ while gen[2] is nd idx
		if gen[2] != nnd
			gen[2] += 1;
		else
			gen[1] += 1;
			gen[2] = 2;
			
			if gen[1] == 2
				ykey = :yᵛ;
				χsY = @view YSOL[ykey].tlvl.nds[1,:];
			else
				ykey = :yⁱ
				χsY = @view YSOL[ykey].tlvl.nds[1,:];
			end
		end

		# Find value of solution directly underneath
		nd = @view EYSOL[ykey].tlvl.nds[:,gen[2]];
		if nd[1] >= χsY[1]
			# Node is overtop the prior t level
			nd0 = γℓvℓ(t₀,nd[1]);
			δτ = nd[2] - nd0[2];

			y₀ = myinterp(χsY,YSOL[ykey].ys,nd[1]);
		else
			# Node is overtop the ∂
			nd0 = [nd[1],0.];
			δτ = nd[2];
			
			if ykey == :yˢ
				y₀ = 0.;
			elseif ykey == :yᵛ
				y₀ = myinterp([ EYSOL[ykey].tlvl.nds[1,1],χsY[1] ],
					      [ ∫λyˢds,YSOL[ykey].ys[1] ],
					      nd[1]);
			elseif ykey == :yⁱ
				y₀ = myinterp([ EYSOL[ykey].tlvl.nds[1,1],χsY[1] ],
					      [ (∫yˢds + ∫Imαyᵛds)*∫βyⁱds,YSOL[ykey].ys[1] ],
					      nd[1]);
			end
		end
		
		# Compute flowfield contribution
		if ykey == :yˢ
			λ₀ = λ(nd0,prm;case=:χτ);
			EYSOL[:yˢ].ys[gen[2]] = y₀ - 1/√(2)*y₀*( λ₀ + ∫βyⁱds )*δτ;
		elseif ykey == :yᵛ
			α₀ = α(nd0,prm;case=:χτ);
			EYSOL[:yᵛ].ys[gen[2]] = y₀ - 1/√(2)*y₀*( α₀ + (1-α₀)*∫βyⁱds )*δτ;
		else
			γ₀ = γ(nd0,prm;case=:χτ);
			EYSOL[:yⁱ].ys[gen[2]] = y₀ - 1/√(2)*γ₀*y₀*δτ;
		end
	end
		
	# Update the intermediate quantities
	βy!(EYSOL[:yⁱ],prm; βy=EYSOL[:βyⁱ]);
	λy!(EYSOL[:yˢ],prm; λy=EYSOL[:λyˢ]);
	Imαy!(EYSOL[:yᵛ],prm; Imαy=EYSOL[:Imαyᵛ]);

	for i=1:nnd
		EYSOL[:λ].ys[i] = λ(EYSOL[:λ].tlvl.nds[:,i],prm;case=:χτ);
		EYSOL[:α].ys[i] = α(EYSOL[:α].tlvl.nds[:,i],prm;case=:χτ);
		EYSOL[:γ].ys[i] = γ(EYSOL[:γ].tlvl.nds[:,i],prm;case=:χτ);
	end

	if flagrt
		return EYSOL
	end

end

# ∂YSOL!
"""
Routine to define the boundary data at the [t=0.] level. 
"""
function ∂YSOL!(prm::DSymFl)
	# Compute dependant params as func of ind
	data!(prm);
	
	# Compute the initial data
	nnd = Int64(prm[:nnd]);

	# Define initial geometry
	tlvl_yˢ = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:Ls],nnd)) );
	tlvl_yᵛ = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:Lv],nnd)) );
	tlvl_yⁱ = Tℓvℓ( 0.,convert(Vector,LinRange(0.,prm[:Li],nnd)) );

	# Define [t=0.] ∂-values
	YSOL = Dict{Symbol,Yℓvℓ}();

	yˢ = Vector{Float64}(undef,nnd);
	yᵛ = zeros(nnd);
	yⁱ = Vector{Float64}(undef,nnd);
	λs = Vector{Float64}(undef,nnd);
	αs = Vector{Float64}(undef,nnd);
	γs = Vector{Float64}(undef,nnd);
	
	# Compute initial values for YSOL
	@inbounds for i=1:nnd
		nd = @view tlvl_yˢ.nds[:,i];
		yˢ[i] = fˢ(nd,prm;case=:χτ);
		λs[i] = λ(nd,prm;case=:χτ);

		nd = @view tlvl_yᵛ.nds[:,i];
		αs[i] = α(nd,prm;case=:χτ);

		nd = @view tlvl_yⁱ.nds[:,i];
		yⁱ[i] = prm[:ρ]*fⁱ(nd,prm;case=:χτ);
		γs[i] = γ(nd,prm;case=:χτ);
	end
	YSOL[:yˢ] = Yℓvℓ(tlvl_yˢ,yˢ);
	YSOL[:yᵛ] = Yℓvℓ(tlvl_yᵛ,yᵛ);
	YSOL[:yⁱ] = Yℓvℓ(tlvl_yⁱ,yⁱ);
	
	YSOL[:λ] = Yℓvℓ(tlvl_yˢ,λs);
	YSOL[:α] = Yℓvℓ(tlvl_yᵛ,αs);
	YSOL[:γ] = Yℓvℓ(tlvl_yⁱ,γs);

	#  βyⁱ,λyˢ,Imαyᵛ
	YSOL[:βyⁱ] = βy!(YSOL[:yⁱ],prm);
	YSOL[:λyˢ] = λy!(YSOL[:yˢ],prm);
	YSOL[:Imαyᵛ] = Imαy!(YSOL[:yᵛ],prm);

	return YSOL
end

# vaxsolver
"""
Integrate the vaccination system by an explicit Euler scheme with adaptive timestep
flagprg says whether to print progress to standard out
"""
function vaxsolver(prm::DSymFl;
		   flagprg::Bool=true)
	# Intialize values and vectors for storing solution 
	updTℓvℓ = [:yˢ,:yᵛ,:yⁱ,:λ,:α,:γ,:βyⁱ,:λyˢ,:Imαyᵛ]; # a list to help reduce mem allocs
	taxis = convert(Vector,0.: prm[:δt] : prm[:T]);
	ntdwn = length(taxis);
	nnd = Int64(prm[:nnd]);
	SOL = Vector{Dict{Symbol,Yℓvℓ}}(undef,ntdwn);
	
	# Setup ∂-[t=0] data
	SOL[1] = ∂YSOL!(prm);

	# Adaptive Euler step requires 4 Yℓvℓ mem locs for writing output
	ynow = deepcopy(SOL[1]); ynext = deepcopy(ynow);
	ymid = deepcopy(ynow); y2xmid = deepcopy(ynow);
	
	pos = 2; # indicates which taxis value is next to surpass
	δt = prm[:δtmax]; # initial guess of adaptive Euler step
	err = Vector{Float64}(undef,nnd); # stores abs nodal errors at [t=tᵢ]
	rnrm = Vector{Float64}(undef,nnd); # stores rel nodal err normalization at [t=tᵢ]
	yaerr = [0.,0.,0.]; # stores the ODE solver absolute error
	yrerr = [0.,0.,0.]; # stores the ODE solver relative error
	yrcap = DSymFl(# stores lower bd caps used to assess rel error
				     # rel err's are wrt to densities over unit length
				     :yˢ=>prm[:rlow]/prm[:Ls],:yᵛ=>prm[:rlow]/prm[:Lv],
				     :yⁱ=>prm[:rlow]/prm[:Li]);
	while pos <= ntdwn
		# Compute the full step
		euler!(δt,ynow,prm;EYSOL=ynext,updTℓvℓ=updTℓvℓ);

		# Compute the two half-steps
		euler!(.5*δt,ynow,prm;EYSOL=ymid,updTℓvℓ=updTℓvℓ);
		euler!(.5*δt,ymid,prm;EYSOL=y2xmid,updTℓvℓ=updTℓvℓ);

		# Compute the yˢ errors
		err[:] = abs.(y2xmid[:yˢ].ys-ynext[:yˢ].ys);
		yaerr[1] = maximum(err);
		mymax!(yrcap[:yˢ],abs.(y2xmid[:yˢ].ys);w=rnrm); 
		yrerr[1] = maximum(err./rnrm);

		# Compute the yᵛ errors
		err[:] = abs.(y2xmid[:yᵛ].ys-ynext[:yᵛ].ys);
		yaerr[2] = maximum(err);
		mymax!(yrcap[:yᵛ],abs.(y2xmid[:yᵛ].ys);w=rnrm); 
		yrerr[2] = maximum(err./rnrm);

		# Compute the yⁱ errors
		err[:] = abs.(y2xmid[:yⁱ].ys-ynext[:yⁱ].ys);
		yaerr[3] = maximum(err);
		mymax!(yrcap[:yⁱ],abs.(y2xmid[:yⁱ].ys);w=rnrm); 
		yrerr[3] = maximum(err./rnrm);

		# Act according to accepting or addapting the t-step
		flagδt = true;
		for i=1:3
			if !flagδt
				break
			end
			# Scale errors by the meas of axis they are densities over
			λnorm = (i==1) ? prm[:Ls] : ( (i==2) ? prm[:Lv] : prm[:Li] )
			flagδt = flagδt&&(
					 (λnorm*yaerr[i]<=prm[:atol])||(yrerr[i]<=prm[:rtol])
					 );
		end

		if !flagδt
			# Error tolerances not attained so adapt and iterate
			δt *= .5;
			
			# Reset the ∫'s computed at midstep so euler! 
			# recomputes at next iteration
			for key in keys(ymid)
				ymid[key].∫yds[:] = [NaN];
			end

			continue
		end
		
		# Error tolerances attained so accept	
		tnow = ynow[:yˢ].tlvl.t₀[1]; tnext = y2xmid[:yˢ].tlvl.t₀[1];
		if tnext >= taxis[pos]	
			# We've passed a downsample t-value and will now store
			posnext = (tnext < taxis[end] ? 
				    myfindfirst(taxis,tnext) : ntdwn + 1 );
			posnext = (pos == posnext ? posnext+1 : posnext);

			# Interpolate the inbetween values
			for i=pos:posnext-1
				SOL[i] = Dict{Symbol,Yℓvℓ}();
				SOL[i][:yˢ] = myinterp([tnow,tnext],[ynow[:yˢ],y2xmid[:yˢ]],taxis[i]);
				SOL[i][:yᵛ] = myinterp([tnow,tnext],[ynow[:yᵛ],y2xmid[:yᵛ]],taxis[i]);
				SOL[i][:yⁱ] = myinterp([tnow,tnext],[ynow[:yⁱ],y2xmid[:yⁱ]],taxis[i]);

				SOL[i][:λ] = myinterp([tnow,tnext],[ynow[:λ],y2xmid[:λ]],taxis[i]);
				SOL[i][:α] = myinterp([tnow,tnext],[ynow[:α],y2xmid[:α]],taxis[i]);
				SOL[i][:γ] = myinterp([tnow,tnext],[ynow[:γ],y2xmid[:γ]],taxis[i]);

				SOL[i][:βyⁱ] = myinterp([tnow,tnext],[ynow[:βyⁱ],y2xmid[:βyⁱ]],taxis[i]);
				SOL[i][:λyˢ] = myinterp([tnow,tnext],[ynow[:λyˢ],y2xmid[:λyˢ]],taxis[i]);
				SOL[i][:Imαyᵛ] = myinterp([tnow,tnext],[ynow[:Imαyᵛ],y2xmid[:Imαyᵛ]],taxis[i]);
			end

			# Define new position
			pos = posnext;

			if flagprg&&(posnext <= ntdwn)
				println("Solving for $posnext/$ntdwn ...");
			end
	
		end
		
		for key in keys(ynow)
			ynow[key].tlvl.t₀[:] = y2xmid[key].tlvl.t₀;
			ynow[key].tlvl.nds[:,:] = y2xmid[key].tlvl.nds;
			ynow[key].tlvl.χrg[:] = y2xmid[key].tlvl.χrg;
			ynow[key].tlvl.τrg[:] = y2xmid[key].tlvl.τrg;
			ynow[key].ys[:] = y2xmid[key].ys;
			ynow[key].∫yds[:] = [NaN];
		end

		for key in keys(ymid)
			ymid[key].∫yds[:] = [NaN];
		end

		# Try a larger time step and continue iteration
		δt = minimum([2*δt,prm[:δtmax]]);
	end

	return taxis,SOL
end

# ∫yʳds
"""
Compute the aggregated recovery time series for the solution output by vaxsolver
Note: prm should be the mutated dictionary output by vaxsolver
"""
function ∫yʳds(taxis::Vector{Float64},SOL::Vector{Dict{Symbol,Yℓvℓ}},prm::DSymFl)
	ntdwn = length(taxis);

	# Initialize starting data
	sol = Vector{Float64}(undef,ntdwn);
	sol[1] = 0.;
	for i=2:ntdwn
		αyᵛ = Yℓvℓ(SOL[i][:yᵛ].tlvl,SOL[i][:yᵛ].ys.*SOL[i][:α].ys);
		γyⁱ = Yℓvℓ(SOL[i][:yⁱ].tlvl,SOL[i][:yⁱ].ys.*SOL[i][:γ].ys);

		flow = ∫line(αyᵛ) + ∫line(γyⁱ);
		sol[i] = prm[:δt]*flow + sol[i-1];
	end

	return sol
end

# savesol
"""
Save the solution into a csv file for plotting and later analysis
"""
function savesol(taxis::Vector{Float64},SOL0::Vector{Dict{Symbol,Yℓvℓ}},
		 prm::DSymFl;
		 fname::String="")
	
	SOL = SOLrefine(prm,SOL0);
	snds = Dict{Symbol,Vector{Float64}}(
		    :yˢ=>SOL[1][:yˢ].tlvl.snds,:yᵛ=>SOL[1][:yᵛ].tlvl.snds,:yⁱ=>SOL[1][:yⁱ].tlvl.snds);
	nnd = length(snds[:yˢ]);
	ntdwn = length(SOL);

	yˢsaxis = reshape(snds[:yˢ],(nnd,1));
	yᵛsaxis = reshape(snds[:yᵛ],(nnd,1));
	yⁱsaxis = reshape(snds[:yⁱ],(nnd,1));

	taxis = reshape(taxis,(ntdwn,1));

	# Write the axis csv's
	CSV.write(fname*"saxis.csv",DataFrame([yˢsaxis yᵛsaxis yⁱsaxis]),header=["ys","yv","yi"]);
	CSV.write(fname*"taxis.csv",DataFrame(taxis),writeheader=false);

	# Write the solution csv's
	Y = Matrix{Float64}(undef,3*nnd,ntdwn);
	
	for i=1:ntdwn
		Y[1:nnd,i] = SOL[i][:yˢ].ys;
		Y[nnd+1:2*nnd,i] = SOL[i][:yᵛ].ys;
		Y[2*nnd+1:3*nnd,i] = SOL[i][:yⁱ].ys;
	end

	CSV.write(fname*"ysol.csv",DataFrame(Y),writeheader=false);

	∫yʳ = ∫yʳds(taxis[:],SOL,prm);
	CSV.write(fname*"yRsol.csv",DataFrame(reshape(∫yʳ,(ntdwn,1))),writeheader=false);
end
