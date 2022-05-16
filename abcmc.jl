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
		 prmrg::DSymVFl=mcdata()[1],
		 prmvary::DSymBool=mcdata()[2],
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
Compute the ℓ² error (squared) between model prediction and data
"""
function ℓerr(ysol::Vector{Solℓvℓ};
	      prm::DSymVFl=data(),
	      yˢ::DStrVFl,yᵛ::DStrVFl,yⁱ::DStrVFl,
	      ram::Vector{Float64}=Vector{Float64}(undef,9))
	valyˢ = 0.0;
	npts = length(yˢ["time"]);
	taxis = [ysol[i].yˢ.tlvl.t₀[1] for i=1:length(ysol)]; 
	@inbounds for k=1:npts
		tnow = yˢ["time"][k];
		ram[1] = yˢ["0-9"][k]; ram[2] = yˢ["10-19"][k]; ram[3] = yˢ["20-29"][k];
		ram[4] = yˢ["30-39"][k]; ram[5] = yˢ["40-49"][k]; ram[6] = yˢ["50-59"][k];
		ram[7] = yˢ["60-69"][k]; ram[8] = yˢ["70-79"][k]; ram[9] = yˢ["80+"][k];
		
		ℓ = myfindfirst(taxis,tnow);
		ℓ = ℓ==1 ? 2 : ℓ;
		ynow = myinterp(tnow,ysol[ℓ-1].yˢ,ysol[ℓ].yˢ);
		ram[1] -= eval(ynow,5.0*365); ram[2] -= eval(ynow,15.0*365); ram[3] -= eval(ynow,25.0*365);
		ram[4] -= eval(ynow,35.0*365); ram[5] -= eval(ynow,45.0*365); ram[6] -= eval(ynow,55.0*365);
		ram[7] -= eval(ynow,65.0*365); ram[8] -= eval(ynow,75.0*365); ram[9] -= eval(ynow,90.0*365);

		# scale density to unit interval
		ram *= prm[:yˢrg][2];

		valyˢ += ram.^2 |> sum;
	end
	# average to one column vector for comparing with yᵛ,yⁱ
	valyˢ *= 1/9;
	
	valyᵛ = 0.0;
	npts = length(yᵛ["time"]);
	@inbounds for k=1:npts
		tnow = yᵛ["time"][k];
		ℓ = myfindfirst(taxis,tnow);
		ℓ = ℓ==1 ? 2 : ℓ;
		ynow = myinterp(tnow,ysol[ℓ-1].yᵛ,ysol[ℓ].yᵛ);

		ram2 = yᵛ["∂"][k]-eval(ynow,0.0);
		# scale density to unit interval
		ram2 *= prm[:yᵛrg][2];

		valyᵛ += ram2^2;

	end

	valyⁱ = 0.0;
	npts = length(yⁱ["time"]);
	@inbounds for k=1:npts
		tnow = yⁱ["time"][k];
		ℓ = myfindfirst(taxis,tnow);
		ℓ = ℓ==1 ? 2 : ℓ;
		ynow = myinterp(tnow,ysol[ℓ-1].yⁱ,ysol[ℓ].yⁱ);

		ram2 = yⁱ["∂"][k]-eval(ynow,0.0);
		# scale density to unit interval
		ram2 *= prm[:yⁱrg][2];

		valyᵛ += ram2^2;
	end

	return valyˢ+valyᵛ+valyⁱ
end

"""
Routine to run abc sampling
"""
function abcrun(nsmp::Int64;
		rng::MersenneTwister=MersenneTwister(),flagrst::Bool=false,
	        δprg::Float64=0.05)
	# Load the data
	dfyˢ = CSV.read("ODH_ys.csv",DataFrame);
	dfyⁱ = CSV.read("ODH_yi.csv",DataFrame);
	dfyᵛ = CSV.read("ODH_yv.csv",DataFrame);

	ODHyˢ=Dict{String,Vector{Float64}}("time"=>(Float64.(dfyˢ[:,:time])),
					"0-9"=>dfyˢ[:,"0-9"],"10-19"=>dfyˢ[:,"10-19"],
					"20-29"=>dfyˢ[:,"20-29"],"30-39"=>dfyˢ[:,"30-39"],
					"40-49"=>dfyˢ[:,"40-49"],"50-59"=>dfyˢ[:,"50-59"],
					"60-69"=>dfyˢ[:,"60-69"],"70-79"=>dfyˢ[:,"70-79"],
					"80+"=>dfyˢ[:,"80+"]);
	ODHyᵛ=Dict{String,Vector{Float64}}("time"=>(Float64.(dfyᵛ[:,:time])),
					"∂"=>dfyᵛ[:,"yv"]);
	ODHyⁱ=Dict{String,Vector{Float64}}("time"=>(Float64.(dfyⁱ[:,:time])),
					"∂"=>dfyⁱ[:,"yi"]);

	# Initialize
	if !flagrst
		prm,vkeys,V = wrtprm();
	else
		dftemp = CSV.read("ABCsmp.csv",DataFrame,header=false);
		vkeys = Symbol.(dftemp[:,1]);
		prm,_ = rdprm(dftemp[!,end],vkeys);
		rng = myloadrng();
	end
	prmrg,prmvary = abcdata();
	S = Matrix{Float64}(undef,length(vkeys),nsmp)

	# Sample
	ram = Vector{Float64}(undef,9);

	println("progress through abc sampling: 0.0/1.0 ...")
	prg = 0.0;
	@inbounds for i=1:nsmp
		abcsmp!(prm;rng=rng,prmrg=prmrg,prmvary=prmvary);
		
		Snow = @view S[:,i];
		wrtprm!(prm,vkeys,Snow);
		CSV.write("ABCsmp.csv",[DataFrame("prm"=>vkeys) DataFrame(S[:,1:i],:auto)],writeheader=false,append=false);
		
		ysol,_ = pdesolve(;prm=prm,flagprg=false); 
		try 
			prm[:ℓerr][1] = ℓerr(ysol;prm=prm,yˢ=ODHyˢ,yᵛ=ODHyᵛ,yⁱ=ODHyⁱ,ram=ram);
		catch
			@warn "simulation failed owing to resolutions and tolerances at sample $i"
		end
		
		wrtprm!(prm,vkeys,Snow);
		CSV.write("ABCsmp.csv",[DataFrame("prm"=>vkeys) DataFrame(S[:,1:i],:auto)],writeheader=false,append=false);
		# Partially save progress and output progess through independent samples and save rng
		while i/nsmp >= prg + δprg
			#CSV.write("ABCsmp.csv",[DataFrame("prm"=>vkeys) DataFrame(S[:,1:i],:auto)],writeheader=false,append=false);
			mysaverng(rng);
			prg += δprg;
			println("progress through abc sampling: $prg/1.0 ...");
		end
	end

	# Save and exit
	println("finished abc sampling ...");
	CSV.write("ABCsmp.csv",[DataFrame("prm"=>vkeys) DataFrame(S,:auto)],writeheader=false,append=false);
	mysaverng(rng);

	return
end
