### A Pluto.jl notebook ###
# v0.12.18

using Markdown
using InteractiveUtils

# ╔═╡ a7b804b0-4923-11eb-37b8-959429a5c0cc
begin
	# hexagonos
	function extInput(x,y,t)
		return 2.0 + 1.0/(0.2^2*pi)*exp(-(x^2+y^2)/0.2^2)
	end

	function kernel(x,y)
		K0 = 0.1
		Kc = 10.0*pi/10
		σ  = 10.0
		ϕ  = [0.0, pi/3.0, 2.0*pi/3.0]
		coskx = 0.0

		for i = 1:3
			coskx += cos(Kc*(x*cos(ϕ[i])+y*sin(ϕ[i])))
		end

		return K0*coskx*exp(-sqrt(x^2+y^2)/σ)
	end

	firingRate(V) = @. 2.0/(1.0+exp(-5.5*(V-3.0)))

	α  = 1.0
	v  = 80.0
	V0 = 2.00083
	L  = 10
	N  = 256
	T  = 1.45
	n  = 145
end

# ╔═╡ 4e1bfe70-54ea-11eb-2365-7758b4d992f8
begin
	using SNFE, Plots
	plotly()
	input = Input2D(α,v,V0,L,N,T,n,extInput,kernel,firingRate);
	prob  = probSNFE(input);
	V     = solveSNFE(prob,[0.5, 0.75, 1.0, 1.25]);
end

# ╔═╡ 605b9ff0-54ea-11eb-39d1-7f90e45e7353
plot(V.x,V.y,V(3),zlims=(2.00083,8))

# ╔═╡ 544ba230-54ec-11eb-0e1a-4d5f9bac4185
begin
	k=[kernel(i,j) for j in V.y, i in V.x]
	plot(V.x,V.y,k)
end

# ╔═╡ Cell order:
# ╠═a7b804b0-4923-11eb-37b8-959429a5c0cc
# ╠═4e1bfe70-54ea-11eb-2365-7758b4d992f8
# ╠═605b9ff0-54ea-11eb-39d1-7f90e45e7353
# ╠═544ba230-54ec-11eb-0e1a-4d5f9bac4185
