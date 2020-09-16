using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings
gr(size=(600,600))

include("../colors.jl")

function beatiful_newton()
  anim = Animation()
  ng = 200
  xg = range(-2, 2, length=ng)
  yg = range(-1.5, 2.5, length=ng)
  NC = 101
  colors = cgrad(color_gradient_list, NC)[1:end]

  f(x) = (x[1] - 1)^2 + 4 * (x[2] - x[1]^2)^2
  f(x,y) = f([x;y])
  ∇f(x) = ForwardDiff.gradient(f, x)
  H(x) = ForwardDiff.hessian(f, x)


  x = [-1.2; 1.0]
  θg = range(0, 2π, length=60)

  ΔN = 10
  state = :fade # :fade => :point

  α = 1.0
  for k = 1:300
    contourf(xg, yg, f, leg=false, grid=false, axis=false, c=colors, levels=100, lw=0)

    fx = f(x)
    g = ∇f(x)
    Hx = H(x)
    m(d) = f(x) + dot(g, d) + dot(d, Hx * d) / 2
    d = -Hx \ g
    y = x + d

    if state == :fade
      Λ, V = eigen(Hx)
      Θ = [cos.(θg)'; sin.(θg)']

      if α > 0
        fz = m(d) * (1 - α) + fx * α
        if fz > fx
          fz = 0.0
        end
        level = 2 * (fz - fx) + dot(-g, d)
        if level > 0
          r = y .+ V' * sqrt.(diagm(1 ./ Λ) * level) * Θ
          plot!(r[1,:], r[2,:], c=:red, lw=3)
          scatter!([x[1]], [x[2]], c=:red, ms=4)
        else
          @info("", level, minimum(eigen(Hx).values))
        end
        α -= 0.05
      else
        α = 1.0
        scatter!([y[1]], [y[2]], c=:red, ms=4)
        state = :point
      end
    elseif state == :point
      plot!([x[1], y[1]], [x[2], y[2]], c=:red, lw=3, l=:arrow)
      scatter!([x[1], y[1]], [x[2], y[2]], c=:red, ms=4)
      if α > 0
        α -= 0.1
      else
        α = 1.0
        state = :fade
        x = y
        if norm(∇f(y)) < 1e-2
          break
        end
      end
    end

    xlims!(xg[1], xg[end])
    ylims!(yg[1], yg[end])
    frame(anim)
  end

  filename = joinpath(@__DIR__, "gif", "beatiful-newton.gif")
  gif(anim, filename, fps=24)
end

beatiful_newton() 