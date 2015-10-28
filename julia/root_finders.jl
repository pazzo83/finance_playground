using Roots
using Polynomials

function incremental_search(f::Function, a::Float64, b::Float64, dx::Float64)
  fa = f(a)
  c = a + dx
  fc = f(c)
  n = 1

  while sign(fa) == sign(fc)
    if a >= b
      return a - dx, n
    end

    a = c
    fa = fc
    c = a + dx
    fc = f(c)
    n += 1
  end

  if fa == 0
    return a, n
  elseif fc == 0
    return c, n
  else
    return (a + c) / 2.0, n
  end
end

function bisection_search(f::Function, a::Float64, b::Float64; tol = 0.1, maxiter = 10)
  n = 1
  while n <= maxiter
    c = (a + b) * 0.5
    if f(c) == 0 || abs(a - b) * 0.5 < tol
      # root found
      return c, n
    end

    n += 1
    if f(c) < 0
      a = c
    else
      b = c
    end
  end

  return c, n
end

function newton_search(f::Function, df::Function, x::Float64; tol = 0.001, maxiter = 100)
  n = 1
  while n <= maxiter
    x1 = x - f(x) / df(x)
    if abs(x1 - x) < tol
      return x1, n
    else
      x = x1
      n += 1
    end
  end

  return nothing, n
end

function secant_search(f::Function, a::Float64, b::Float64; tol = 0.001, maxiter = 100)
  n = 1

  while n <= maxiter
    c = b - f(b) * ((b - a) / (f(b) - f(a)))
    if abs(c - b) < tol
      return c, n
    end

    a = b
    b = c
    n += 1
  end

  return nothing, n
end

# y(x) = x^3 + 2.0x^2 - 5.0
# dy(x) = 3.0x^2 + 4.0x
# ddy(x) = 6.0x + 4.0
#
# root, iterations = incremental_search(y, -5.0, 5.0, 0.001)
# root, iterations = bisection_search(y, -5.0, 5.0; tol = 0.00001)
# root, iterations = newton_search(y, dy, 5.0; tol = 0.00001)
# root, iterations = secant_search(y, -5.0, 5.0; tol = 0.00001)
#
# # using the Roots module
# # fzero
# fzero(y, 5.0) # one guess
# fzero(y, [-5.0, 5.0]) # two guesses
# fzeros(y) # no guesses
#
# #newton
# newton(y, 5.0) # no deriv
# newton(y, dy, 5.0) # deriv
#
# #halley
# halley(y, 5.0) # no deriv
# halley(y, dy, ddy, 5.0) # two derivs
#
# #secant
# secant_method(y, -5.0, 5.0)
#
# #using Polynomials
# poly_y = Poly([-5.0, 0, 2.0, 1.0])
# roots(poly_y)
