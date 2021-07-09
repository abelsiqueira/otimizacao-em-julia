### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 0dbc216a-7b23-4d44-8f3c-57a19345a2c4
using PlutoUI, Plots, PyCall, Colors, Distributions

# ╔═╡ 56574e89-fe54-403d-af9f-79bf6aa17c83
begin
	acoes = ["VALE3.SA","PETR4.SA","ITSA4.SA","MGLU3.SA","AZUL4.SA"] |> sort
	pd = pyimport("pandas")
	yf = pyimport("yfinance")
	data = yf.download(;tickers=acoes, period="1y", group_by="column")[:Close].values
end

# ╔═╡ e93606aa-25cd-4c4b-8b80-82afa90f57ff
begin
	colors = Colors.distinguishable_colors(length(acoes) + 2)[3:end]
	colors = reshape(colors, 1, :)
	labels = reshape(acoes, 1, :)
	plt_args = (:leg => :topleft, :lab => labels, :c => colors)
end

# ╔═╡ 1b8da58f-875d-4379-bed6-88dd4d10211c
plot(data; plt_args...)

# ╔═╡ 60659d95-805d-4ac2-80b1-dc36688c9ece
plot(data ./ data[1,:]'; plt_args...)

# ╔═╡ df02a8e7-1e3c-4c3c-aa90-d29d9fa864bb
retornos = (data[2:end,:] - data[1:end-1,:]) ./ data[1:end-1,:]

# ╔═╡ e74b5c14-8f41-49aa-ac4d-0e298c743e12
plot(retornos; plt_args...)

# ╔═╡ a08ea546-4c9b-4bb9-86c2-80f901b16b5e
histogram(retornos, bins=50, opacity=0.7; plt_args...)

# ╔═╡ 2b3f6965-201f-44e7-907c-9f51ac55a021
begin
	μM = reshape(mean(retornos, dims=1), :, 1)
	σM = cov(retornos)
	histogram(retornos, bins=50, normalize=true, opacity=0.5; plt_args...)
	for i = 1:length(acoes)
		μᵢ = μM[i]
		σᵢ = sqrt(σM[i,i])
		d = Normal(μᵢ, σᵢ)
		plot!(x -> pdf(d, x), μᵢ - 3σᵢ, μᵢ + 3σᵢ, lw=3; c=colors[i], lab="")
	end
	plot!()
end

# ╔═╡ 347e92e6-9b00-41d9-872f-f6683961cdf1
md"""
x1 = $(@bind x1 Slider(0:0.05:1, show_value=true))
x2 = $(@bind x2 Slider(0:0.05:1, show_value=true))
x3 = $(@bind x3 Slider(0:0.05:1, show_value=true))
x4 = $(@bind x4 Slider(0:0.05:1, show_value=true))
"""

# ╔═╡ 2e7038a2-504d-4267-bf48-2daa20a6ef57
begin
	using LinearAlgebra
	
	x = [x1; x2; x3; x4]
	x = [x; 1 - sum(x)]
	if x[end] ≥ 0
		cart_ret = retornos * x
		cart_val = [1.0; cumprod(1 .+ cart_ret)]

		μ = mean(cart_ret)
		σ = std(cart_ret)

		plt = plot(layout=grid(2, 2), size=(800, 600), leg=false)
		for i = 1:length(acoes)
			μᵢ = μM[i]
			σᵢ = sqrt(σM[i,i])
			d = Normal(μᵢ, σᵢ)
			plot!(plt[1,1], x -> pdf(d, x), μᵢ - 3σᵢ, μᵢ + 3σᵢ, lw=3; c=colors[i], lab=labels[i])
			plot!(plt[1,1], [μᵢ, μᵢ], [0, pdf(d, μᵢ)], c=colors[i], lab="")
			plot!(plt[1,1], [μᵢ - σᵢ, μᵢ + σᵢ], pdf(d, μᵢ - σᵢ) * [1, 1], c=colors[i], lab="")
		end
		normal = Normal(μ, σ)
		plot!(plt[1,1], x -> pdf(normal, x), μ - 3σ, μ + 3σ, lw=3, c=:black, lab="cart")
		plot!(plt[1,1], [μ, μ], [0, pdf(normal, μ)], c=:black, lab="")
		plot!(plt[1,1], [μ - σ, μ + σ], pdf(normal, μ - σ) * [1, 1], c=:black, lab="")

		scatter!(plt[1,2], μM', sqrt.(diag(σM)'), xlabel="retorno", ylabel="risco"; plt_args...)
		scatter!(plt[1,2], [μ], [σ], c=:black, lab="", leg=false)
		xlims!(plt[1,2], 0, maximum(μM))
		ylims!(plt[1,2], 0, maximum(sqrt.(diag(σM))))
		# μθ = []
		# σθ = []
		# for θ = range(0, 1, length=50)
		# 	xθ = [θ; 1 - θ]
		# 	push!(μθ, dot(μM, xθ))
		# 	push!(σθ, sqrt(xθ' * σM * xθ))
		# end
		# plot!(plt[1,2], μθ, σθ, c=:black, lab="")

		plot!(plt[2,1], data ./ data[1,:]'; plt_args...)
		plot!(plt[2,1], cart_val, c=:black, lab="", lw=3)

		plot!(plt[2,2], retornos; plt_args...)
		plot!(plt[2,2], cart_ret, c=:black, lab="")

		plt
	end
end

# ╔═╡ Cell order:
# ╠═0dbc216a-7b23-4d44-8f3c-57a19345a2c4
# ╠═56574e89-fe54-403d-af9f-79bf6aa17c83
# ╠═e93606aa-25cd-4c4b-8b80-82afa90f57ff
# ╠═1b8da58f-875d-4379-bed6-88dd4d10211c
# ╠═60659d95-805d-4ac2-80b1-dc36688c9ece
# ╠═df02a8e7-1e3c-4c3c-aa90-d29d9fa864bb
# ╠═e74b5c14-8f41-49aa-ac4d-0e298c743e12
# ╠═a08ea546-4c9b-4bb9-86c2-80f901b16b5e
# ╠═2b3f6965-201f-44e7-907c-9f51ac55a021
# ╟─347e92e6-9b00-41d9-872f-f6683961cdf1
# ╠═2e7038a2-504d-4267-bf48-2daa20a6ef57
