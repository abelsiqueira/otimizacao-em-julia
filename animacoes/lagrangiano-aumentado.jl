using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,600))

include("../colors.jl")

function lagrangiano_aumentado()
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
  b = [2; 3; 2; 1; 1; 1]
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
  # f(x,y) = (x - 1)^2 + 2 * (y - x^2)^2
  f(x) = f(x...)
  
  x = [1.4; 1.4]
  X = reshape(x, 2, 1)
  y = zeros(length(b))
  
  for ρ = 2 .^ range(-8, 6, length=60)
    ϕ(x) = f(x) + dot(y, A * x - b) + ρ * sum(max.(0, A * x - b).^2)
    ϕ(x,y) = ϕ([x;y])
    
    ∇f(x) = ForwardDiff.gradient(f, x)
    ∇ϕ(x) = ForwardDiff.gradient(ϕ, x)
    Hϕ(x) = ForwardDiff.hessian(ϕ, x)
    iter = 0
    while norm(∇f(x) - A' * y) > 1e-1
      Δx = -Hϕ(x) \ ∇ϕ(x)
      iter += 1
      x += Δx
      y += max.(A * x - b, 0) * ρ
      X = [X x]
      iter > 3 && break
    end

    plot(leg=false, grid=false, axis=false)
    contourf!(xg, yg, ϕ, c=colors, lw=0, levels=NC)

    scatter!(F[1,:], F[2,:], c=:black, ms=5)
    plot!(F[1,:], F[2,:],    c=:black, lw=4)
    plot!(F[1,:], F[2,:],    c=:white, lw=2)
    scatter!(F[1,:], F[2,:], c=:white, ms=3)

    plot!(X[1,:], X[2,:], c=:red, lw=2, l=:arrow)

    xlims!(max(xg[1], x[1] - 1 / ρ^2), min(xg[end], x[1] + 1 / ρ^2))
    ylims!(max(yg[1], x[2] - 1 / ρ^2), min(yg[end], x[2] + 1 / ρ^2))

    frame(anim)
  end
  filename = joinpath(@__DIR__, "gif", "lagrangiano-aumentado.gif")
  gif(anim, filename, fps=15)
end

lagrangiano_aumentado()