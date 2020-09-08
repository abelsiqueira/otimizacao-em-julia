using Colors, Plots, LinearAlgebra, LaTeXStrings
pyplot(size=(600,600))

include("../colors.jl")

function quadraticas_metodos()
  anim = Animation()
  ng = 100
  xg = range(-2, 2, length=ng)
  yg = range(-2, 2, length=ng)
  NC = 21
  colors = cgrad(color_gradient_list, NC)
  V = [2.0 1; -1 2] / √5
  A = V' * diagm(0 => [6.0; 0.5]) * V
  f(x) = dot(x, A * x)
  f(x,y) = f([x;y])

  for θ = range(0, 2π, length=24*4)[1:end-1]
    θ += π / 2
    x = 1.5 * [cos(θ); sin(θ)]
    plot(leg=false, grid=false, axis=false)
    contourf!(xg, yg, f, leg=false, levels=100, c=colors)
    scatter!([x[1]], [x[2]], c=:white, ms=6)

    # Conjugate gradient
    r = d = -A * x
    α = dot(d, d) / dot(d, A * d)
    y = x + α * d
    rt = -A * y
    β = dot(rt, rt) / dot(r, r)
    v = rt + β * d
    v = dot(rt, rt) / dot(v, A * v) * v
    z = y + v
    plot!([x[1], y[1], z[1]], [x[2], y[2], z[2]], c=:red, lw=2)

    # Cauchy
    X = zeros(2, 50)
    X[:,1] = x
    for k = 1:size(X,2)-1
      d = -A * X[:,k]
      d = dot(d, d) / dot(d, A * d) * d
      X[:,k+1] = X[:,k] + d
    end
    plot!(X[1,:], X[2,:], c=:black, lw=2, l=:dash)

    xlims!(xg[1], xg[end])
    ylims!(yg[1], yg[end])
    frame(anim)
  end

  filename = joinpath(@__DIR__, "gif", "quadraticas-metodos.gif")
  gif(anim, filename, fps=24)
end

quadraticas_metodos()