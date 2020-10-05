using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,600))

include("../colors.jl")

function project(x, c, ∇c, Hc)
  y = copy(x)
  λ = 0.0
  iter = 0
  while norm(y - x + λ * ∇c(y)) + norm(λ * c(y)) + norm(max(c(y), 0)) > 1e-8
    K = [I + λ * Hc(y)  ∇c(y); ∇c(y)' 0.0]
    rhs = [y - x + λ * ∇c(y); c(y)]
    Δ = -K \ rhs
    y += Δ[1:2]
    λ += Δ[3]
    iter += 1
    if iter > 10
      break
    end
  end
  return y
end


function restricoes_convexo()
  anim = Animation()
  NC = 100

  xg = range(-4, 4, length=NC)
  yg = range(-4, 4, length=NC)

  colors = cgrad(color_gradient_list, NC)[1:end]

  α = 1.5
  # c(x,y) = (9α^2 * x^2 + 18α^2 * y^2 + 3α^3 * x * y^2 + α^2 * y^2 - 144) / 100
  # c(x,y) = 4 * x^2 + y^2 - 3
  c(x,y) = log(exp(-x-2) + exp(-2y-2) + exp(3x + 4y - 4) + exp(-6x - 5y - 2)) - 1.5
  c(x) = c(x...)

  ∇c(x) = ForwardDiff.gradient(c, x)
  Hc(x) = ForwardDiff.hessian(c, x)

  for θ in range(0, 2π, length=120)
  # for θ in [9π / 12]
    # x = [cos(θ); sin(θ)] * 3.8
    γ = [cos(θ); sin(θ)] * 3
    x = zeros(2)
    f(x,y) = cos(θ)^2 * (x - γ[1])^2 + sin(θ)^2 * (y - γ[2])^2 + cos(θ) * sin(θ) * (x - γ[1]) * (y - γ[2])
    f(x) = f(x...)
    ∇f(x) = ForwardDiff.gradient(f, x)
    Hf(x) = ForwardDiff.hessian(f, x)

    plot(leg=false, grid=false, axis=false)
    contourf!(xg, yg, f, levels=NC, c=colors, lw=0, leg=false)
    contour!(xg, yg, c, levels=[0.0], lw=2, c=:black, leg=false)

    iter = 0
    g = ∇f(x)
    α = dot(g, g) / dot(g, Hf(x) * g)
    y = project(x - α * g, c, ∇c, Hc)
    while norm(y - x) > 1e-2
      if x - α * g != y
        plot!(x[1] .+ [0, -α * g[1]], x[2] .+ [0, -α * g[2]], c=:orange, l=:arrow)
        plot!([x[1] - α * g[1], y[1]], [x[2] - α * g[2], y[2]], c=:green, l=:arrow)
      end
      plot!([x[1], y[1]], [x[2], y[2]], c=:white, m=(3,stroke(1,:black)), l=:arrow)
      x = y
      g = ∇f(x)
      α = dot(g, g) / dot(g, Hf(x) * g)
      y = project(x - α * g, c, ∇c, Hc)
      iter += 1
      if iter > 10
        break
      end
    end

    # L(x) = c(y) + dot(∇c(y), y - x)
    # L(x,y) = L([x;y])

    # contour!(xg, yg, L, levels=[0.0], c=:red, l=:dash)

    frame(anim)
  end
  filename = joinpath(@__DIR__, "gif", "restricoes-convexos.mp4")
  mp4(anim, filename, fps=15)
  filename = joinpath(@__DIR__, "gif", "restricoes-convexos.gif")
  filenametrim = joinpath(@__DIR__, "gif", "restricoes-convexos-trim.gif")
  gif(anim, filename, fps=15)
  run(`convert $filename -trim +repage $filenametrim`)
end

restricoes_convexo()