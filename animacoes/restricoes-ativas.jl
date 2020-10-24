using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,600))

include("../colors.jl")

function restricoes_ativas()
  anim = Animation()
  NC = 200

  xg = range(-1.2, 1.5, length=NC)
  yg = range(-1.2, 1.5, length=NC)

  colors = cgrad(color_gradient_list, NC)[1:end]

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

  # for t = -1.2:0.1:1.2
  for θ = range(0, 2π, length=180)
    γ = [cos(θ); sin(θ)] * 1.4
    # f(x,y) = cos(θ)^2 * (x - γ[1])^2 + sin(θ)^2 * (y - γ[2])^2 + cos(θ) * sin(θ) * (x - γ[1]) * (y - γ[2])
    # f(x,y) = (x - γ[1])^2 + (y - γ[2])^2
    f(x,y) = 5.0 * (x - 0.5 - 0.8 * sin(θ))^2 + (y - 0.8 - 0.5 * cos(θ))^2 * (5 + 4 * sin(θ))
    # f(x,y) = 0.3 * (x - cos(θ))^2 + 0.3 * 4 * (y - sin(θ) * x^2)^2 + (x-1)^2 + (y-1)^2
    f(x) = f(x...)
    ∇f(x) = ForwardDiff.gradient(f, x)
    Hf(x) = ForwardDiff.hessian(f, x)

    x = [-1.0; -1.0]
    X = reshape(x, 2, 1)
    W = findall(abs.(A * x - b) .< 1e-6)

    gx = ∇f(x)
    y = zeros(m)
    y[W] = -A[W,:]' \ gx
    solved = norm(gx + A' * y) < 1e-6 && all(y .≥ 0)
    iter = 0
    tired = iter > 100

    while !(solved || tired)
      if norm(gx + A' * y) < 1e-6
        # Candidate
        if any(y .< 0)
          i = findfirst(y .< 0)
          setdiff!(W, i)
          y = zeros(m)
          y[W] = -A[W,:]' \ gx
        end
      else
        K = [Hf(x) A[W,:]'; A[W,:] zeros(length(W), length(W))]
        rhs = [gx + A[W,:]' * y[W]; A[W,:] * x - b[W]]
        Δx = (-K \ rhs)[1:2]
        α = 1.0
        wall = 0
        for i = setdiff(1:m, W)
          β = A[i,:]' * Δx
          if β > 0
            αt = (b[i] - A[i,:]' * x) / β
            if αt < α
              α = αt
              wall = i
            end
          end
        end
        x += α * Δx
        X = [X x]
        if wall != 0
          push!(W, wall)
        end
      end
      gx = ∇f(x)
      y = zeros(m)
      y[W] = -A[W,:]' \ gx

      solved = norm(gx + A' * y) < 1e-6 && all(y .≥ 0)
      iter += 1
      tired = iter > 100
    end

    plot(leg=false, grid=false, axis=false)
    contourf!(xg, yg, f, c=colors, lw=0, levels=NC)

    scatter!(F[1,:], F[2,:], c=:black, ms=5)
    plot!(F[1,:], F[2,:], c=:black, lw=4)
    plot!(F[1,:], F[2,:], c=:white, lw=2)
    scatter!(F[1,:], F[2,:], c=:white, ms=3)

    plot!(X[1,:], X[2,:], c=:red, lw=2, l=:arrow)

    frame(anim)
  end
  filename = joinpath(@__DIR__, "gif", "restricoes-ativas.gif")
  gif(anim, filename, fps=15)
end

restricoes_ativas()