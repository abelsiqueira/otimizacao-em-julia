using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,600))

include("../colors.jl")

function pontos_interiores()
  anim = Animation()
  NC = 200

  xg = range(-1.2, 1.5, length=NC)
  yg = range(-1.2, 1.5, length=NC)

  colors = cgrad(color_gradient_list, NC)[1:end]
  darker_colors = weighted_color_mean.(0.25, RGBA(0.5, 0.5, 0.5, 1.0), colors)

  Random.seed!(0)
  # A = round.(Int, [10 * randn(m, 2) ./ rand(m); -I])
  # b = round.(Int, [10 * randn(m) ./ rand(m); 0; 0])
  A = [3 0.2;
  2 6;
  -2.5 1;
  1 -1;
  -1 0;
  0 -1]
  b = [4; 8; 2; 1; 1; 1]
  m = length(b)
  # c(x) = sum(max.(A * x - b, 0))

  feasible_set = Vector{Float64}[]
  for i = 1:m-1, j = i+1:m
    x = A[[i,j],:] \ b[[i,j]]
    if all(A * x .≤ b .+ 1e-6)
      push!(feasible_set, x)
    end
  end
  xc = sum(feasible_set) / length(feasible_set)
  θ = [atan(f - xc ...) for f in feasible_set]
  J = sortperm(θ)
  F = hcat([feasible_set[J]; [feasible_set[J[1]]]]...)

  f(x,y) = (x - 0.5)^2 + (y - 1.3)^2 * 5
  f(x) = f(x...)
  
  x = [-0.99; -0.99]
  X = reshape(x, 2, 1)
  
  for μ = 2 .^ range(12, -12, length=60)
    ϕ(x) = all(A * x .≤ b) ? f(x) - μ * sum(log.(b - A * x)) : NaN
    ϕ(x,y) = ϕ([x;y])
    
    ∇ϕ(x) = ForwardDiff.gradient(ϕ, x)
    Hϕ(x) = ForwardDiff.hessian(ϕ, x)
    iter = 0
    while norm(∇ϕ(x)) > max(μ * 1e-1, 1e-1)
      Δx = -Hϕ(x) \ ∇ϕ(x)
      α = 1.0
      while any(A * (x + α * Δx) .> b)
        α /= 2
      end
      iter += 1
      x += α * Δx
      X = [X x]
      iter > 10 && break
    end

    plot(leg=false, grid=false, axis=false)
    contourf!(xg, yg, f, c=darker_colors, lw=0, levels=range(0.01, 40, length=250))
    contourf!(xg, yg, ϕ, c=colors, lw=0, levels=NC)

    scatter!(F[1,:], F[2,:], c=:black, ms=5)
    plot!(F[1,:], F[2,:], c=:black, lw=4)
    plot!(F[1,:], F[2,:], c=:white, lw=2)
    scatter!(F[1,:], F[2,:], c=:white, ms=3)

    plot!(X[1,:], X[2,:], c=:red, lw=2, l=:arrow)

    frame(anim)
  end
  filename = joinpath(@__DIR__, "gif", "pontos-interiores.gif")
  gif(anim, filename, fps=15)
end

pontos_interiores()