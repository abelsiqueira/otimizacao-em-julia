using Colors, ForwardDiff, NLPModels, NLPModelsIpopt, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,600))

include("../colors.jl")

function movement_in_field()
  anim = Animation()
  NC = 20

  xg = range(-1.5, 2.5, length=NC)
  yg = range(-1.5, 2.5, length=NC)

  colors = cgrad(color_gradient_list, NC)[1:end]

  P(x,y) = -y
  Q(x,y) = x
  M = maximum(norm([P(xi,yi); Q(xi,yi)]) for xi=xg, yi=yg)
  Δ = sqrt((xg[2] - xg[1])^2 + (yg[2] - yg[1])^2)

  N = 1000
  # model = Model()
  # @variable(model, xg[1] ≤ x[0:N] ≤ xg[end])
  # @variable(model, yg[1] ≤ y[0:N] ≤ yg[end])
  # JuMP.register(model, :P, 2, P, autodiff=true)
  # JuMP.register(model, :Q, 2, Q, autodiff=true)
  # JuMP.fix(x[0], 0.0, force=true)
  # JuMP.fix(y[0], 0.0, force=true)
  # JuMP.fix(x[N], 1.0, force=true)
  # JuMP.fix(y[N], 1.0, force=true)
  f(x,y) = sum(P(x[n],y[n]) * (x[n+1] - x[n]) + Q(x[n],y[n]) * (y[n+1] - y[n]) for n = 1:N) +
           1.0 * sum(sqrt((x[n+1] - x[n])^2 + (y[n+1] - y[n])^2) for n = 1:N)
  f(x) = f(x[1:N+1], x[N+2:end])
  ℓ = [0.0; fill(xg[1], N-1); 1.0; 0.0; fill(yg[1], N-1); 1.0]
  u = [0.0; fill(xg[end], N-1); 1.0; 0.0; fill(yg[end], N-1); 1.0]
  nlp = ADNLPModel(f, [range(0, 1, length=N+1); range(0, 1, length=N+1)], ℓ, u)
  # nlp = MathOptNLPModel(model)

  for k = 1:1

    plot(leg=false, grid=false, axis=false, ticks=false)
    for xi = xg, yi = yg
      v = [P(xi,yi); Q(xi,yi)]
      v = v / M * Δ
      plot!(xi .+ [0, v[1]], yi .+ [0, v[2]], l=:arrow, c=:black)
    end

    x = nlp.meta.x0
    fx = obj(nlp, x)
    iter = 0
    g = grad(nlp, x)
    π = norm(x - clamp.(x - g, ℓ, u))
    while !(π < 1e-6 || iter > 5000)
      if norm(g) > 0.1
        g *= 0.1 / norm(g)
      end
      α = 1.0
      xt = clamp.(x - α * g, ℓ, u)
      ft = obj(nlp, xt)
      while !(ft < fx + 1e-2 * dot(g, xt - x))
        α /= 2
        if α < 1e-6
          break
        end
        xt = clamp.(x - α * g, ℓ, u)
        ft = obj(nlp, xt)
      end
      x = xt
      fx = ft
      g = grad(nlp, x)
      π = norm(x - clamp.(x - g, ℓ, u))
      iter += 1
      # plot!(x[1:N+1], x[N+2:end], c=:blue, lw=1)
    end
    println("π = $π, iter = $iter")
    plot!(x[1:N+1], x[N+2:end], c=:red, lw=2, m=(2))

    xlims!(xg[1], xg[end])
    ylims!(yg[1], yg[end])

    frame(anim)
  end

  filename = joinpath(@__DIR__, "gif", "movement-in-field.gif")
  gif(anim, filename, fps=15)
end

movement_in_field()