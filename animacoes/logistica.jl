using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,400))

include("../colors.jl")

function logistica()
  anim = Animation()
  ng = 200
  β1g = range(-2, 2, length=ng)
  β2g = range(0, 6, length=ng)
  NC = 101
  colors = cgrad(color_gradient_list, NC)[2:end]

  n = 50
  Random.seed!(0)
  x = sort(2 * rand(n) .- 1)
  y = 1 ./ (1 .+ exp.(-x + randn(n) * 0.45)) .> 0.5
  h(β, x) = 1 / (1 + exp(-β[1] - β[2] * x))
  L(β) = -sum(
    y[i] * log(h(β, x[i])) + (1 - y[i]) * log(1 - h(β, x[i])) for i = 1:n
  )

  N = 50
  state = :vertical
  θ = 0.0
  while state != :none
    layout = @layout [
      a{0.2h}
      grid(1, 2)
    ]
    p = plot(
      layout=layout,
      grid=false,
      axis=false,
      leg=false
    )
    annotate!(p[1], [0.0], [-1], text(L"L(\beta) = -\sum_i y_i \ln\ \left(\!h(\beta, x_i)\right) + (1 - y_i) \ln\ (1 - h(\beta, x_i))", 15))
    annotate!(p[1], [0.0], [1], text(L"h(\beta, x) = (1 + e^{-\beta_0 - \beta_1 x})^{-1}", 15))
    xlims!(p[1], -1.2, 1.2)
    ylims!(p[1], -1.2, 1.2)
    scatter!(p[2], x, y, leg=false, yticks=[0,1])
    contourf!(p[3], β1g, β2g, (a,b) -> L([a;b]), levels=200, leg=false, c=colors, lw=0)

    β = [0.0; 3.0]
    if state == :vertical
      β[2] += sin(θ) * 2.5
      θ += 2π / N
      if θ > 2π
        state = :horizontal
        θ = 0.0
      end
    elseif state == :horizontal
      β[1] += sin(θ) * 1.5
      θ += 2π / N
      if θ > 2π
        state = :none
        θ = 0.0
      end
    end
    plot!(p[2], x -> h(β, x), extrema(x)..., c=:red, lw=2)
    scatter!(p[3], [β[1]], [β[2]], c=:red, ms=4)

    xlabel!(p[2], "x")
    ylabel!(p[2], "y")
    xlabel!(p[3], "β₀")
    ylabel!(p[3], "β₁")

    frame(anim)
  end

  filename = joinpath(@__DIR__, "gif", "logistica.gif")
  gif(anim, filename, fps=24)
end

logistica()