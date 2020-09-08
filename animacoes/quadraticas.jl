using Colors, Plots, LinearAlgebra, LaTeXStrings
gr()

include("../colors.jl")

function conicas()
  anim = Animation()
  ng = 40
  xg = range(-2, 2, length=ng)
  yg = range(-2, 2, length=ng)
  NC = 21
  colors = cgrad(color_gradient_list, NC)
  colors_c = cgrad(repeat(color_gradient_list, 8), 8NC)
  zℓ = -8

  θt = range(0, 2π, length=ng)
  N = 40
  ΔN = 40
  for α in [
    range(-1.0, 0.0, length=ΔN);
    range(0.0, 1.0, length=ΔN);
    range(1.0, -1.0, length=ΔN)
  ]
    layout = @layout [
      grid(1,3){0.15h}
      grid(1,3){0.6h}
      grid(1,3){0.15h}
      a{0.1h}
    ]
    p = plot(layout=layout, size=(600,600), leg=false, axis=false, grid=false)
    F = [
      (x,y) -> α * x^2 + y^2
      (x,y) -> α * x^2 - y^2
      (x,y) -> x^2 + α * y
    ]

    for (i,t) in enumerate([L"\alpha x^2 + y^2", L"\alpha x^2 - y^2", L"x^2 + \alpha y"])
      annotate!(p[i], [0], [0], text(t, 25))
      xlims!(p[i], -1, 1)
      ylims!(p[i], -1, 1)
    end

    for k = 4:6
      surface!(p[k], xg, yg, F[k-3], c=colors, camera=(20,20), framestyle=:none, grid=false)
      # contour!(p[k], xg, yg, (x,y) -> F[k-3](x,y) / 8 - 7, c=colors_c, levels=20)
      xlims!(p[k], extrema(xg)...)
      ylims!(p[k], extrema(yg)...)
      zlims!(p[k], zℓ, 8)
    end

    if α > 0
      annotate!(p[7], [0], [0.5], text("min", 20))
      annotate!(p[7], [0], [-0.5], text("único", 20))
    elseif α == 0
      annotate!(p[7], [0], [0.5], text("min", 20))
      annotate!(p[7], [0], [-0.5], text(L"\infty", 20))
    else
      annotate!(p[7], [0], [0.5], text("sela", 20))
      annotate!(p[7], [0], [-0.5], text("único", 20))
    end
    xlims!(p[7], -1, 1)
    ylims!(p[7], -1, 1)

    if α > 0
      annotate!(p[8], [0], [0.5], text("sela", 20))
      annotate!(p[8], [0], [-0.5], text("único", 20))
    elseif α == 0
      annotate!(p[8], [0], [0.5], text("max", 20))
      annotate!(p[8], [0], [-0.5], text(L"\infty", 20))
    else
      annotate!(p[8], [0], [0.5], text("max", 20))
      annotate!(p[8], [0], [-0.5], text("único", 20))
    end
    xlims!(p[8], -1, 1)
    ylims!(p[8], -1, 1)

    if α != 0
      annotate!(p[9], [0], [0.0], text(L"\nexists", 20))
    else
      annotate!(p[9], [0], [0.0], text(L"\infty", 20))
    end
    xlims!(p[9], -1, 1)
    ylims!(p[9], -1, 1)

    plot!(p[10], [-1, 1], [0, 0], c=:black, lw=2)
    for x = -1:1
      plot!(p[10], [x, x], [-1e-1, 1e-1], c=:black, lw=2)
      annotate!(p[10], [x], [-2e-1], text("$x"))
    end
    scatter!(p[10], [α], [0.0], c=:red, ms=3)
    scatter!(p[10], [α], [0.1], c=:red, ms=6, m=:diamond)
    ylims!(p[10], -0.3, 0.2)

    frame(anim)
    if α * (α^2 - 1) == 0
      for i = 1:N
        frame(anim)
      end
    end
  end

  filename = joinpath(@__DIR__, "gif", "quadraticas.gif")
  gif(anim, filename, fps=24)
end

conicas()