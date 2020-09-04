using Colors, Plots, LinearAlgebra, JuMP, Cbc
gr(size=(600,600))

include("../colors.jl")

function alicebob3d()
  a = 10.0
  L = [30; 20]
  H = [4; 3]
  Q = 3.0
  U = [0.5; 0.3]

  model = Model(Cbc.Optimizer)
  @variable(model, x[1:2] ≥ 0, Int)
  @variable(model, y ≥ 0, Int)
  @objective(model, Max, dot(L, x) - a * y)
  @constraint(model, dot(H, x) ≤ 80)
  @constraint(model, dot(U, x) .≤ Q + y)
  optimize!(model)

  x = value.(x)
  y = value(y)
  @info("",
    x, y, dot(L, x) - a * y
  )

  x1max = 22.0
  x2max = 30.0
  ymax = 20.0

  anim = Animation()
  N = 240
  colors = cgrad(color_gradient_list, 200)
  for (i,c) = enumerate(range(1, 800, length=N))
    c = min(c, 530)
    color = colors[ceil(Int, c / 530 * length(colors))]
    plot(leg=false)

    for k = 1:2:c
      plot!(
        [k / L[1], 0],
        [0, k / L[2]],
        [0, 0],
        c=colors[ceil(Int, k / 530 * length(colors))]
      )
    end

    # axis
    plot!([0, x1max], [0, 0], [0, 0], c=:gray, lw=1)
    plot!([0, 0], [0, x2max], [0, 0], c=:gray, lw=1)
    plot!([0, 0], [0, 0], [0, ymax],  c=:gray, lw=1)
    # H[1] * x[1] + H[2] * x[2] ≤ 80
    plot3d!(
      [80 / H[1], 80 / H[1], 0,         0,         80 / H[1]],
      [0,         0,         80 / H[2], 80 / H[2], 0],
      [0,         ymax,      ymax,      0,         0],
      c=:gray, l=:dash
    )
    plot3d!([80 / H[1], 0], [0, 80 / H[2]], [ymax, ymax], c=:gray, lw=1)
    # U[1] * x[1] + U[2] * x[2] - y ≤ Q
    plot3d!(
      [Q / U[1], 0,        0,                 (Q + ymax) / U[1], Q / U[1]],
      [0,        Q / U[2], (Q + ymax) / U[2], 0,                 0],
      [0,        0,        ymax,              ymax,              0],
      c=:gray, l=:dash
    )

    # Intersecção
    plot3d!(
      [0,                    80 / H[1]],
      [80 / H[2],            0],
      [U[2] * 80 / H[2] - Q, U[1] * 80 / H[1] - Q],
      c=:black, lw=2
    )

    plot3d!(
      [Q / U[1], 80 / H[1], 80 / H[1], 0, 0, 0, Q / U[1]],
      [0, 0, 0, 80 / H[2], 80 / H[2], Q / U[2], 0],
      [0, U[1] * 80 / H[1] - Q, ymax, ymax, U[2] * 80 / H[2] - Q, 0, 0],
      c=:black, lw=2
    )


    # L[1] * x[1] + L[2] * x[2] - a * y = c
    plot3d!(
      [c / L[1], 0,        0,                     (c + a * ymax) / L[1], c / L[1]],
      [0,        c / L[2], (c + a * ymax) / L[2], 0,                     0],
      [0,        0,        ymax,                  ymax,                  0],
      c=color, lw=2,
    )
    # Intersecção com U[1] * x[1] + U[2] * x[2] - y = Q
    x1 = (c - a * Q) / (L[1] - a * U[1])
    x2 = (c - a * Q) / (L[2] - a * U[2])
    plot3d!(
      [0,             x1],
      [x2,            0],
      [U[2] * x2 - Q, U[1] * x1 - Q],
      c=color, lw=2
    )
    # Intersecção com H[1] * x[1] + H[2] * x[2] = 80
    x1 = 80 / H[1]
    x2 = 80 / H[2]
    plot3d!(
      [0,             x1],
      [x2,            0],
      [(L[2] * x2 - c) / a, (L[1] * x1 - c) / a],
      c=color, lw=2
    )

    xlims!(0, x1max)
    ylims!(0, x2max)
    zlims!(0, ymax)
    frame(anim)
  end
  filename = joinpath(@__DIR__, "gif", "alicebob3d.jl")
  gif(anim, filename, fps=30)
end

alicebob3d()