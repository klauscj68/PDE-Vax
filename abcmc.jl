using Random,Distributions

"""
Output parameter ranges and bool flags to say which parameters are being 
varied in abc sampling. Note these ranges are before application of data!
to fix parameters that depend on others.
"""
function abcdata()
	prmrg = Dict{Symbol,Vector{Float64}}();
	prmvary = Dict{Symbol,Bool}();
	
	# α parameters
	prmrg[:αL] = [7.0,21.0];
	prmvary[:αL] = false;

	prmrg[:αeff] = [0.7,1.0];
	prmvary[:αeff] = true;

	# β parameters
	#  θ-scale param not included because it is determined by the f^i
	#  initial condition
	prmrg[:βα]=[5.0,20.0];
	prmvary[:βα] = true;

	# γ parameters
	prmrg[:γθ]=[1.0,18.0];
	prmvary[:γθ] = true;

	prmrg[:γα]=[2.0,8.0];
	prmvary[:γα] = true;

	# λ parameters
	prmrg[:λ]=[0.1,10.0];
	prmvary[:λ] = false;

	# mass of t₀ infected
	prmrg[:ρ] = [0.001,0.05];
	prmvary[:ρ] = true;

	return prmrg,prmvary
end

"""
Sample varied parameters like priors
"""
function abcsmp!(prm::DSymVFl;
		 prmrg::DSymVFl=abcdata()[1],
		 prmvary::DSymBool=abcdata()[2],
		 rng::MersenneTwister=MersenneTwister())
	flagfd = false;

	while !flagfd
		for key in keys(prmvary)
			if prmvary[key]
				prm[key][1] = prmrg[key][1]+rand(rng)*(prmrg[key][2]-prmrg[key][1]);
			end
		end

		# Enforce dependent relationships
		data!(prm);

		# Check if recovery is less than 3 weeks
		γ = Weibull(prm[:γα][1],prm[:γθ][1]);
		γqt = quantile(γ,0.999);
		if γqt > 21.0
			continue
		end

		# Check if infectious transmission is shorter than recovery
		β = Weibull(prm[:βα][1],prm[:βθ][1]);
		γμ = mean(γ); βμ = mean(β);
		if βμ > γμ
			continue
		end
		flagfd = true;
	end
	
end

""" 
Pregenerate abcsamples in specified batch sizes to run embarassingly 
parallel on supercomputer
Each batch will have nsmp's each for a total of nsmp*nbatch many samples
"""
function abcpregen(nsmp::Int64,nbatch::Int64;
		   prm::DSymVFl=data(),
		   prmrg::DSymVFl=abcdata()[1],
		   prmvary::DSymBool=abcdata()[2],
		   rng::MersenneTwister=MersenneTwister())
	vkeys = [mykey for mykey in keys(prm)];
	Stemp = Matrix{Float64}(undef,length(vkeys),nsmp);

	gen = [1,0]
	@inbounds for k=1:nsmp*nbatch
		if gen[2]!=nsmp
			gen[2]+=1
		else
			gen[1]+=1; gen[2]=1;
		end

		abcsmp!(prm;prmrg=prmrg,prmvary=prmvary,rng=rng);
		Snow = @view Stemp[:,gen[2]];
		wrtprm!(prm,vkeys,Snow)

		if gen[2]==nsmp
			dftemp = DataFrame(hcat(vkeys,Stemp),:auto);
			id=gen[1];
			CSV.write("ABCsmp$id.csv",dftemp,writeheader=false);
		end
	end 
	mysaverng(rng);
end
"""
Compute the ℓ² error (squared) between model prediction and data
daily incidence
"""
function ℓerr(ysol::Vector{Solℓvℓ};
	      prm::DSymVFl=data(),
	      df::DataFrame=CSV.read("ODH_snipdaily.csv",DataFrame))
	npts = nrow(df);
	taxis = [ysol[i].yˢ.tlvl.t₀[1] for i=1:length(ysol)]; 	

	val = 0.0;
	# Compute effective population size for extracting model pred
	# daily incidence
	#  Total infections during this period
	kT = sum(df[!,"daily_confirm"]);

	#  ∫yⁱdt at s=0 by trapezoidal rule
	∫yⁱdt = 0.0;
	@inbounds for k=2:length(taxis)
		∫yⁱdt += (ysol[k].yⁱ.ys[1]+ysol[k-1].yⁱ.ys[1])/2*(taxis[k]-taxis[k-1]);
	end
	neff = kT/∫yⁱdt;

	# Compute difference in daily incidence between model prediction and observed
	@inbounds for k=1:npts
		tnow = k-1.0;
		ℓ = myfindfirst(taxis,tnow);
		ℓ = ℓ==1 ? 2 : ℓ;
		ynow = myinterp(tnow,ysol[ℓ-1].yⁱ,ysol[ℓ].yⁱ);

		val += ( df[!,"daily_confirm"][k]-neff*eval(ynow,0.0) )^2;
	end

	return √(val/npts)
end

"""
Routine to run abc sampling. fsmp is the file output by abcpregen
that you are using for samples. pos is optional flag to say where to
pick up at (for ex if sheet had been partially run through)
"""
function abcrun(fsmp::String;
	        δprg::Float64=0.05,
		fODH::String="ODH_snipdaily.csv",
		pos::Int64=1)
	# Load the samples
	dfS = CSV.read(fsmp,DataFrame,header=false);
	nsmp = ncol(dfS)-1;
	vkeys = Symbol.(dfS[:,1]);

	# find pos of :ℓerr within vkeys
	posℓ = 1;
	while vkeys[posℓ]!=:ℓerr
		posℓ+=1;
	end

	# Load the data
	dfODH = CSV.read(fODH,DataFrame);
	select!(dfODH,["time","daily_confirm"]);

	println("progress through abc sampling: 0.0/1.0 ...")
	prg = 0.0;
	@inbounds for i=(pos+1):ncol(dfS)	
		prm = rdprm(dfS[!,i],vkeys)[1];
		ysol,_ = pdesolve(;prm=prm,flagprg=false);

		try 	
			dfS[posℓ,i] = ℓerr(ysol;prm=prm,df=dfODH);
		catch
			dfS[posℓ,i] = NaN;
			@warn "simulation failed owing to resolutions and tolerances at sample $i"
		end
		
		CSV.write(fsmp,dfS,writeheader=false,append=false);
		# Partially save progress and output progess
		while i/nsmp >= prg + δprg
			#CSV.write(fsmp,dfS,writeheader=false,append=false);	
			prg += δprg;
			println("progress through abc sampling: $prg/1.0 ...");
		end
	end

	# Save and exit
	println("finished abc sampling ...");
	CSV.write(fsmp,dfS,writeheader=false,append=false);

	return
end
