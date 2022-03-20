# Collection of routines for computing with 1d S^r_d(Δ) bsplines on an interval
using Plots,Measures
gr();
VFl = Vector{Float64};
VecVw = Union{Vector{Float64},
	      SubArray{Float64, 1, Vector{Float64}, Tuple{UnitRange{Int64}}, true},
	     };

"""
Define a class for an extended knot sequence for S^r_d(Δ)
"""
struct Knots
	Δ::VFl
	Δext::VFl
	r::Int64
	d::Int64
	dim::Int64

	function Knots(Δ::VecVw,r::Int64,d::Int64)
		nknt = length(Δ);
		nintknt = nknt-2;	
		# Compute extended knot partition
		Δext=Vector{Float64}(undef,2*(d+1)+(d-r)*nintknt);
		v = @view Δext[1:d+1]; fill!(v,Δ[1]);
		v = @view Δext[end-d:end]; fill!(v,Δ[end]);
		@inbounds for i=2:nknt-1
			v = @view Δext[d+1+(i-2)*(d-r)+1:d+1+(i-1)*(d-r)]; fill!(v,Δ[i]);
		end

		# Compute the dimension of S^r_d(Δ)
		dim = (nintknt+1)*(d+1)-nintknt*(r+1);

		return new(Δ[:],Δext,r,d,dim)
	end

	function Knots(Δ::VecVw,Δext::VecVw,r::Int64,d::Int64)
		nknt = length(Δ);
		nintknt = nknt-2;

		# Compute the dimension
		dim = (nintknt+1)*(d+1)-nintknt*(r+1);

		return new(Δ[:],Δext[:],r,d,dim)
	end
end
function Base.show(io::IO,kn::Knots)
	for fld in fieldnames(Knots)
		val = getfield(kn,fld);
		flds = string(fld);
		println(io,flds*": $val");
	end
end

"""
Define a class 1d bsplines on S^r_d(Δ)
"""
struct Bspl1d
	c::VFl
	knots::Knots

	function Bspl1d(c::VecVw,kn::Knots)
		if length(c)<kn.dim
			cnew = vcat(c,fill(0.0,kn.dim-length(c)));
		elseif length(c)>kn.dim
			cnew = c[1:kn.dim];
		else
			cnew = c[:];
		end

		return new(cnew,kn)
	end
	function Bspl1d(c::VecVw,Δ::VecVw,r::Int64,d::Int64)
		kn = Knots(Δ,r,d);

		return Bspl1d(c,kn)
	end
end
function Base.show(io::IO,bspl::Bspl1d)
	val = bspl.c;
	println("c: $val");
	println(bspl.knots);
end

"""
auxilliary routine needed to evaluate a spline
Find the first proper interval [,) of extended knot sequence containing
the point. Returns l=dim of S^r_d if point equals right end
"""
function findinterval(knots::Knots,t::Float64)
	@assert (t>=knots.Δ[1])&&(t<=knots.Δ[end]) "query point must be in interval";
	ℓ = 1;
	while (ℓ<=knots.dim-1)&&(knots.Δext[ℓ+1]<=t)
		ℓ += 1;
	end

	return ℓ
end

"""
Evaluate a spline at the given query point by De'Castlejeau
"""
function eval!(bspl1d::Bspl1d,t::Float64;
	      ram::VecVw=Vector{Float64}(undef,bspl1d.knots.d+1))
	npts = bspl1d.knots.d+1;

	# initialize first row, note loop cycles backwards over row
	ℓ = findinterval(bspl1d.knots,t);
	@inbounds for k=1:npts
		ram[npts+1-k] = bspl1d.c[ℓ+1-k];
	end

	# Run De'Castlejeau
	#  In Schumaker, bottom index is ℓ+1-j and his formula is for c^[j+1]_i
	#  while you below are working with c^[i-1]_(ℓ+1-j)
	gen = [npts-1,1];
	while true
		if (gen[2] == gen[1]+1)&&(gen[1] > 1)
			gen[1] -= 1; gen[2] = 1;
		elseif (gen[2] == gen[1]+1)&&(gen[1] == 1)
			break
		end
		i = npts+1-gen[1]; j = gen[2];
		y0 = bspl1d.knots.Δext[ℓ+1-j]; y = bspl1d.knots.Δext[ℓ+1-j+npts-(i-1)];
		ram[end+1-j] = y == y0 ? 0.0 : (t-y0)/(y-y0)*ram[end+1-j] + (y-t)/(y-y0)*ram[end-j];

		gen[2] += 1;
	end

	return ram[end]
end
function eval!(bspl1d::Bspl1d,vt::VecVw;
	       ramDC::VecVw=Vector{Float64}(undef,bspl1d.knots.d+1),
	       ramval::VecVw=Vector{Float64}(undef,length(vt)))
	@inbounds for i=1:length(vt)
		ramval[i] = eval!(bspl1d,vt[i];ram=ramDC);
	end
end

"""
Routine to compute the derivative of a spline
"""
function ∂(bspl1d::Bspl1d;
	   ram::VecVw=Vector{Float64}(undef,bspl1d.knots.dim-1))
	Δext = bspl1d.knots.Δext[2:end-1];
	@inbounds for i=1:bspl1d.knots.dim-1
		ram[i] = (bspl1d.knots.Δext[i+bspl1d.knots.d+1]==bspl1d.knots.Δext[i+1] ? 
			    0.0 : bspl1d.knots.d*(bspl1d.c[i+1]-bspl1d.c[i])/(bspl1d.knots.Δext[i+bspl1d.knots.d+1]-bspl1d.knots.Δext[i+1]));
	end
	return Bspl1d(ram,Knots(bspl1d.knots.Δ,Δext,bspl1d.knots.r-1,bspl1d.knots.d-1))
end

"""
Routine to compute the antiderivative of a spline. Note that the spline returned
is for the definite integral from the left end of knot sequence to future evaluation
point.
"""
function ∫(bspl1d::Bspl1d;
	   ram::VecVw=Vector{Float64}(undef,bspl1d.knots.dim+1))
	Δext = [bspl1d.knots.Δ[1];bspl1d.knots.Δext[:];bspl1d.knots.Δ[end]];
	ram[1] = 0.0;
	gen = [2,1];
	while gen[1]<=bspl1d.knots.dim+1
		if gen[2] == 1
			ram[gen[1]] = 0.0
		end
		i = gen[1]; j = gen[2];
		ram[i] += bspl1d.c[j]*(bspl1d.knots.Δext[bspl1d.knots.d+j+1]-bspl1d.knots.Δext[j]);

		if gen[2]<gen[1]-1
			gen[2] += 1;
		else
			gen[1] += 1; gen[2] = 1;
		end
	end
	ram *= 1/(bspl1d.knots.d+1);
	return Bspl1d(ram,Knots(bspl1d.knots.Δ,Δext,bspl1d.knots.r+1,bspl1d.knots.d+1))
end

"""
Output all S^r_d(Δ) basis splines in an ordered vector
"""
function basis(knots::Knots;
	          ram::VecVw=fill(0.0,knots.dim))
	n = knots.dim;
	V = Vector{Bspl1d}(undef,n);

	@inbounds for i=1:n
		ram[i] = 1.0;
		V[i] = Bspl1d(ram,knots);
		ram[i] = 0.0;
	end

	return V
end

"""
Compute the least squares approximating spline in S^r_d(Δ) for the data
"""
function lsqspl(x::VFl,y::VFl,knots::Knots;
		ramG::Matrix{Float64}=Matrix{Float64}(undef,knots.dim,knots.dim),
		ramf::VecVw=Vector{Float64}(undef,knots.dim),
		ramDC::VecVw=Vector{Float64}(undef,knots.d+1),
		basis::Vector{Bspl1d}=basis(knots))
	npts = length(y);

	# Build Gram matrix
	gen = [1,1,0];
	@inbounds for k=1:knots.dim*knots.dim*npts
		if gen[3]<npts
			gen[3]+=1;
		elseif gen[2]<knots.dim
			gen[2]+=1; gen[3]=1;
		else
			gen[1]+=1; gen[2]=1; gen[3]=1;
		end
		i₁=gen[1];i₂=gen[2];i₃=gen[3];
		if i₃==1
			ramG[i₁,i₂] = 0.0;
		end
		
		ℓ=findinterval(knots,x[i₃]);
		if (i₁<ℓ-knots.d)||(i₁>ℓ)||(i₂<ℓ-knots.d)||(i₂>ℓ)
			continue
		end

		ramG[i₁,i₂] += eval!(basis[i₁],x[i₃];ram=ramDC)*eval!(basis[i₂],x[i₃];ram=ramDC);

	end

	# Build forcing term
	gen = [1,0];
	@inbounds for k=1:knots.dim*npts
		if gen[2]<npts
			gen[2] += 1;
		else
			gen[1] += 1; gen[2] = 1;
		end
		i=gen[1]; j=gen[2];
		
		if j==1
			ramf[i]=0.0;
		end

		ramf[i] += y[j]*eval!(basis[i],x[j];ram=ramDC);
	end

	# Compute interpolating spline
	cnew = ramG\ramf;

	return Bspl1d(cnew,knots)
end

"""
Simple arithmetic operations on splines. Standard addition and subtraction
should be on splines in same S^r_d(Δ) space
"""
function Base.:+(bspl1d1::Bspl1d,bspl1d2::Bspl1d)
	cnew = bspl1d1.c+bspl1d2.c;

	return Bspl1d(cnew,bspl1d1.knots)
end
function Base.:+(bspl1d::Bspl1d,k::Real)
	cnew = bspl1d.c .+ k;

	return Bspl1d(cnew,bspl1d.knots)
end
function Base.:+(k::Real,bspl1d::Bspl1d)
	return bspl1d+k
end
function Base.:-(bspl1d1::Bspl1d,bspl1d2::Bspl1d)
	cnew = bspl1d1.c-bspl1d2.c;

	return Bspl1d(cnew,bspl1d1.knots)
end
function Base.:-(bspl1d::Bspl1d,k::Real)
	return bspl1d+(-k)
end
function Base.:-(k::Real,bspl1d::Bspl1d)
	return k+Bspl1d(-1*bspl1d.c,bspl1d.knots)
end
function Base.:*(λ::Real,bspl1d::Bspl1d)
	cnew = λ*bspl1d.c;

	return Bspl1d(cnew,bspl1d.knots)
end
function Base.:*(bspl1d::Bspl1d,λ::Real)
	return λ*bspl1d
end


"""
Plotting routine for splines, including a single spline, vector of splines, and
the basis splines
"""
function RecipesBase.plot(bspl1d::Bspl1d;npts=1001)
	xval = LinRange(bspl1d.knots.Δ[1],bspl1d.knots.Δ[end],npts);
	yval = [eval!(bspl1d,x) for x in xval];

	plot(xval,yval,xlabel="x",ylabel="y",linewidth=3,margin=5mm);
	plot!()
end
function RecipesBase.plot(bspl1ds::Vector{Bspl1d};npts=1001)
	n = length(bspl1ds);
	p1 = plot();
	@inbounds for i=1:n
		xval = LinRange(bspl1ds[i].knots.Δ[1],bspl1ds[i].knots.Δ[end],npts);
		yval = [eval!(bspl1ds[i],x) for x in xval];
		
		plot!(p1,xval,yval,xlabel="x",ylabel="y",linewidth=3,margin=5mm);
	end
	plot!()
end
function RecipesBase.plot(knots::Knots;npts=1001)
	# Create vector of basis splines
	v = Vector{Bspl1d}(undef,knots.dim);

	c = fill(0.0,knots.dim);
	@inbounds for i=1:knots.dim
		c[i] = 1.0;
		v[i] = Bspl1d(c,knots);
		c[i] = 0.0;
	end

	plot(v;npts=npts)
	plot!()
end
function RecipesBase.plot(Δ::VecVw,r::Int64,d::Int64;npts=1001)
	knots = Knots(Δ,r,d);

	plot(knots;npts=1001)
	plot!()
end
