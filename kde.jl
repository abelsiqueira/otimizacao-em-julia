### A Pluto.jl notebook ###
# v0.15.0

using Markdown
using InteractiveUtils

# ╔═╡ 3c772464-d919-11eb-3bea-65c621b7dd97
using Dates, Plots, Distributions, Statistics, MarketData, PlutoUI

# ╔═╡ a0036e97-5a74-495c-afed-2030b0d08a1c
begin
	opt = YahooOpt(
		period1 = DateTime("2019-01-01"),
		period2 = DateTime("2020-12-31"),
		interval = "1d",
	)
	data = yahoo("AZUL4.SA", opt).Close
	data = values(data) |> skipmissing |> collect
end

# ╔═╡ 04681200-409e-497c-bdb9-dabc12b45d86
begin
	# retornos = (data[2:end] - data[1:end-1]) ./ data[1:end-1]
	retornos = rand(Normal(0, 0.02), 200)
	retornos[1:10:end] .= 0.1 .+ randn(length(1:10:200)) * 0.005
	retornos[2:10:end] .= -0.1 .+ randn(length(1:10:200)) * 0.005
end

# ╔═╡ 8c62b9d3-541a-40a8-b7aa-57c961786eea
begin
	n = length(retornos)
	μ = mean(retornos)
	σ = std(retornos)
	D = Normal(μ, σ)
	h = 1.06 * σ / n^(1/5)
	K(x) = exp(-x^2 / 2) / √(2π)
	Kh(x) = K(x / h) / h
	f(x) = sum(Kh(x - xi) for xi in retornos) / n
	plot(grid=false, axis=false, leg=false, xticks=[], yticks=[], background_color=RGBA(0, 0, 0, 0), size=(800,600))
	histogram!(retornos, bins=50, normalize=true)
	plot!(f, c=:blue, lw=2)
	plot!(x -> pdf(D, x), c=:red, lw=2)
	# title!("h = $h")
	png("kde-1d")
	plot!()
end

# ╔═╡ 299d54bf-5b3b-4d4c-8adc-b619e3d81618
begin
	using LinearAlgebra, Random
	let
		Random.seed!(0)
		n = 400
		X = [randn(div(n, 4), 2) .* [0.2; 0.1]';
			 randn(div(n, 4), 2) .* [0.2; 0.1]' .+ [2; 1]';
			 randn(div(n, 4), 2) * 0.2 .+ [1; 0]';
			 randn(div(n, 4), 2) .* [0.2; 0.6]' .+ [0.5; 0.5]']
		
		μ = [0.0; 0.0]
		Σ = [1.0 0.0;
			 0.0 1.0]
		D = MultivariateNormal(μ, Σ)
		K(x) = pdf(D, x)
		h = 0.25
		Kh(x) = K(x / h) / h
		H = [0.05 0.01; 0.01 0.02]
		G = cholesky(H).L
		KH(x) = K(G \ x) / √det(H)
		f(x) = sum(KH(x - X[i,:]) for i = 1:n) / n
		
		plot(grid=false, axis=false, leg=false, xticks=[], yticks=[], background_color=RGBA(0, 0, 0, 0), size=(800,600))
		scatter!(X[:,1], X[:,2], m=(3,:lightblue,stroke(1,:gray)), leg=false, axis_ratio=:equal)
		contour!(
			range(minimum(X[:,1])-0.1, maximum(X[:,1])+0.1, length=100),
			range(minimum(X[:,2])-0.1, maximum(X[:,2])+0.1, length=100),
			(x, y) -> f([x; y]),
			lw=1
		)
		png("kde-contour-2d")
		plot!()
	end
end

# ╔═╡ bd29df35-35dc-4253-8217-9c25ae14297b
begin
	normal = Normal(0, 1)
	IK(x) = cdf(normal, x)
	IKh(x) = IK(x / h)
	F(x) = sum(IKh(x - xi) for xi in retornos) / n
	plot(F, extrema(retornos)..., leg=false)
	plot!(x -> cdf(D, x), extrema(retornos)..., c=:red)
end

# ╔═╡ b659c7d2-c002-488f-9fdb-00933b9a314f
begin
	let
		Random.seed!(0)
		n = 400
		X = [randn(div(n, 4), 2) .* [0.2; 0.1]';
			 randn(div(n, 4), 2) .* [0.2; 0.1]' .+ [2; 1]';
			 randn(div(n, 4), 2) * 0.2 .+ [1; 0]';
			 randn(div(n, 4), 2) .* [0.2; 0.6]' .+ [0.5; 0.5]']
		
		μ = [0.0; 0.0]
		Σ = [1.0 0.0;
			 0.0 1.0]
		D = MultivariateNormal(μ, Σ)
		K(x) = pdf(D, x)
		h = 0.25
		Kh(x) = K(x / h) / h
		H = [0.05 0.01; 0.01 0.02]
		G = cholesky(H).L
		KH(x) = K(G \ x) / √det(H)
		f(x) = sum(KH(x - X[i,:]) for i = 1:n) / n
		
		plot(grid=true, axis=false, leg=false, xticks=[], yticks=[], zticks=[], background_color=RGBA(0, 0, 0, 0), size=(800,600), camera=(50, 50))
		scatter!(X[:,1], X[:,2], zeros(n), m=(3,:lightblue,stroke(1,:gray)), leg=false, axis_ratio=:equal)
		surface!(
			range(extrema(X[:,1])..., length=100),
			range(extrema(X[:,2])..., length=100),
			(x, y) -> f([x; y]),
			opacity=0.5,
			lw=1
		)
		png("kde-2d-surface")
		plot!()
	end
end

# ╔═╡ Cell order:
# ╠═3c772464-d919-11eb-3bea-65c621b7dd97
# ╠═a0036e97-5a74-495c-afed-2030b0d08a1c
# ╠═04681200-409e-497c-bdb9-dabc12b45d86
# ╠═8c62b9d3-541a-40a8-b7aa-57c961786eea
# ╟─bd29df35-35dc-4253-8217-9c25ae14297b
# ╠═299d54bf-5b3b-4d4c-8adc-b619e3d81618
# ╠═b659c7d2-c002-488f-9fdb-00933b9a314f
