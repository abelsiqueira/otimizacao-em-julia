using Colors, ForwardDiff, Plots, LinearAlgebra, LaTeXStrings, Random
gr(size=(600,600))

include("../colors.jl")

function exactsearch(f)
  a = 0.0
  b = 1.0
  γ = (3 - sqrt(5)) / 2
  while f(b) > f(a) - 1e-4 * (b - a)^2
    b = a + (b - a) * γ
  end
  c = a + (b - a) * γ

  while f(c) ≥ f(b)
    a, c = c, b
    b = a + (c - a) / γ
    if b > 1000
      error("ilimitado")
    end
  end

  d = b - (b - a) * γ

  while abs(b - a) > 1e-4
    if f(c) < f(d)
      a, b, d = a, d, c
      c = a + (b - a) * γ
    else
      a, b, c = c, b, d
      d = b - (b - a) * γ
    end
  end

  if f(c) < f(d)
    return c
  else
    return d
  end
end

function trustregion(B, g, Δ)
  λℓ = 0.0
  λ = 1.0
  d = -B \ g
  if norm(d) ≤ Δ
    return d
  end
  d = -(B + λ * I) \ g
  while norm(d) > Δ
    λℓ, λ = λ, 10λ
    d = -(B + λ * I) \ g
  end
  λu = λ
  λ = (λℓ + λu) / 2
  d = -(B + λ * I) \ g
  while abs(norm(d) - Δ) > 1e-4
    if norm(d) > Δ
      λℓ = λ
    else
      λu = λ
    end
    λ = (λℓ + λu) / 2
    d = -(B + λ * I) \ g
  end
  if norm(d) > Δ
    d = d * Δ / norm(d)
  end
  return d
end

function beatiful_quasi_newton()
  ng = 200
  xg = range(-2, 2, length=ng)
  yg = range(-2, 4, length=ng)
  NC = 101
  colors = cgrad(color_gradient_list, NC)[1:end]

  # f(x) = sqrt( (x[1] - 1)^2 + 4 * (x[2] - x[1]^2)^2 ) / 5 + 0.01 *  sin(dot(x,x)) * dot(x,x) + sum((x .- 1).^4) * 1e-3
  Random.seed!(5)
  A = randn(5, 2) ./ rand(5)
  f(x) = log(sum(exp.(A * (x .- 1))))
  f(x,y) = f([x;y])
  ∇f(x) = ForwardDiff.gradient(f, x)
  H(x) = ForwardDiff.hessian(f, x)

  anim = Animation()

  # for method in [:Newton], control in [:tr]
  for method in [:Gradiente, :Newton, :BFGS, :DFP], control in [:ls, :tr]
    x = [-1.5; -1.0]
    B = if method == :Newton
      H(x)
    else
      Matrix(1.0I, 2, 2)
    end
    X = reshape(x, 2, 1)
    g = ∇f(x)
    Δ = norm(g)
    iter = 0
    while norm(g) > 1e-8 && iter < 10
      local d
      updated = false
      if control == :ls
        iter += 1
        updated = true
        d = -B \ g
        if dot(d, g) > 0
          d = -d
        end
        # α = exactsearch(α -> f(x + α * d))
        α = 1.0
        while f(x + α * d) ≥ f(x) + α * 1e-2 * dot(d, g)
          α /= 2
          if α < 1e-12
            break
          end
        end
        d = α * d
        X = [X  x + d]
        x = x + d
        g = ∇f(x)
      else
        d = trustregion(B, g, Δ)
        Ared = f(x) - f(x + d)
        Pred = -dot(g, d) - dot(d, B * d) / 2
        ρ = Ared / Pred
        if dot(g, d) ≥ 0
          λ = eigen(B).values
          @info("TR fail", dot(g, d), Ared, Pred, ρ, λ[1], λ[2])
          break
        end
        if ρ > 1e-2
          iter += 1
          X = [X  x + d]
          x = x + d
          g = ∇f(x)
          if ρ > 0.75
            Δ *= 4
          end
          updated = true
        else
          Δ /= 4
        end
      end
      if method == :Newton
        B = H(x)
      elseif method == :BFGS
        y = ∇f(x + d) - g
        if dot(y, d) > 0 && updated
          B = B - (B * d) * (B * d)' / dot(d, B * d) + y * y' / dot(y, d)
        end
      elseif method == :DFP
        y = ∇f(x + d) - g
        ρ = 1 / dot(y, d)
        if dot(y, d) > 0 && updated
          B = (I - ρ * y * d') * B * (I - ρ * d * y') + ρ * y * y'
        end
      end
      @info("$iter", g, norm(g))
    end

    cols = size(X, 2) - 1
    for j = 1:cols, α in range(0, 1, length=5)
    # for j = size(X, 2)-1, α = 1.0
      Xq = X[:,1:j+1]
      Xs = [X[:,1:j]  X[:,j] * (1 - α) + X[:,j+1] * α]
      # Everything up to j
      contourf(xg, yg, f, c=colors, levels=100, lw=0, leg=false, axis=false, grid=false)
      plot!(Xq[1,:], Xq[2,:], m=(4,:black,stroke(0)), l=:arrow, c=:black, lw=3)
      plot!(Xs[1,:], Xs[2,:], m=(3,:black,stroke(1,:white)), l=:arrow, c=:white, lw=2)
      title = "$method com " * (control == :ls ? "busca linear" : "região de confiança")
      # title!(title)
      λ = (j + α - 1) / cols
      # annotate!([xg[1] + xg[end]]/2, [3.3],
      #           text(title, weighted_color_mean(λ, colorant"black", colorant"white"), 16))
      annotate!([xg[1] + xg[end] + 0.2 * (1 - 2λ)]/2, [3.2], text(title, :black, 16))
      annotate!([xg[1] + xg[end]]/2, [3.3], text(title, :white, 18))
      frame(anim)
    end
    # for k = 1:12
    #   frame(anim)
    # end
  end

  filename = joinpath(@__DIR__, "gif", "beatiful-quasi-newton.gif")
  gif(anim, filename, fps=24)
  # png(filename)
end

beatiful_quasi_newton()