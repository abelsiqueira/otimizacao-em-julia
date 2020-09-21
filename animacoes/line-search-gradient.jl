using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings
gr(size=(600,600))

include("../colors.jl")

function line_search_gradient()
  anim = Animation()
  ng = 200
  xg = range(-2, 2, length=ng)
  yg = range(-1.5, 2.5, length=ng)
  NC = 101
  colors = cgrad(color_gradient_list, NC)[1:end]

  f(x) = ( sqrt((x[1] - 1)^2 + 4 * (x[2] - x[1]^2)^2) + (x[1]^2 + x[2]^2) * sin(x[1]^2 + x[2]^2) * 0.2 ) / 3
  f(x,y) = f([x;y])
  ∇f(x) = ForwardDiff.gradient(f, x)

  x = [-1.2; 1.0]

  ΔN = 15
  state = :tangent # :tangent => :move

  X = zeros(2,0)
  X = [X x]
  α = 0.0
  Δα = 0.5
  αmax = 3.0
  for k = 1:400
    contourf(xg, yg, f, leg=false, grid=false, axis=false, c=colors, levels=100, lw=0)

    fx = f(x)
    g = ∇f(x)
    n = [g[2]; -g[1]] / norm(g)
    many = false

    if size(X, 1) > 0
      plot!(X[1,:], X[2,:], c=:white, lw=2, m=(5,stroke(0),:white))
    end

    if state == :tangent
      contour!(xg, yg, f, c=:red, levels=[fx], lw=2)
      scatter!([x[1]], [x[2]], c=:white, ms=5)
      if α > 0
        plot!(
          [x[1] - n[1] * α, x[1] + n[1] * α],
          [x[2] - n[2] * α, x[2] + n[2] * α],
          c=:white, lw=2, l=:dash
        )
        plot!(
          [x[1], x[1] - α / αmax * g[1]],
          [x[2], x[2] - α / αmax * g[2]],
          c=:white,
          l=:arrow,
          lw=2
        )
      end
      α += Δα
      if α > αmax
        state = :move
        α = 1.0
        many = true
      end
    elseif state == :move
      plot!(
        [x[1] - n[1] * αmax, x[1] + n[1] * αmax],
        [x[2] - n[2] * αmax, x[2] + n[2] * αmax],
        c=:white, lw=2, l=:dash
      )
      plot!(
        [x[1], x[1] - α * g[1]],
        [x[2], x[2] - α * g[2]],
        c=:white,
        l=:arrow, lw=2
      )
      contour!(xg, yg, f, c=:red, levels=[fx], lw=2)
      contour!(xg, yg, f, c=:yellow, levels=[f(x - α * g)], lw=2)
      if f(x - α * g) > fx - 0.1 * α * dot(g, g)
        α *= 0.99
      else
        many = true
        state = :tangent
        x = x - α * g
        X = [X x]
        α = 0.0
      end
    end

    xlims!(xg[1], xg[end])
    ylims!(yg[1], yg[end])
    frame(anim)
    if many
      for _ in 1:ΔN
        frame(anim)
      end
    end
  end

  filename = joinpath(@__DIR__, "gif", "line-search-gradient.gif")
  gif(anim, filename, fps=24)
end

line_search_gradient() 