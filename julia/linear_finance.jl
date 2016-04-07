# for CAPM
calc_beta(stockreturns::Vector{Float64}, mktreturns::Vector{Float64}) = linreg(stockreturns, mktreturns)

# for Arbitrage Pricing Theory model
ols_estimate(x_vals::Matrix{Float64}, y_vals::Vector{Float64}) = linreg(x_vals, y_vals)

function linalg_solve()
  # 2a + b + c = 4
  # a + 3b + 2c = 5
  # a = 6
  A = [2 1 1; 1 3 2; 1 0 0]
  B = [4, 5, 6]

  # Ax = B -> x = A^-1 * B
  return A \ B
end

function linalg_solve_lu()
  # 2a + b + c = 4
  # a + 3b + 2c = 5
  # a = 6
  A = [2 1 1; 1 3 2; 1 0 0]
  B = [4, 5, 6]
  LU = lufact(B)

  return LU \ B
end

function linalg_solve_ch()
  A = [10 -1 2 0; -1 11 -1 3; 2 -1 10 -1; 0 3 -1 8]
  B = [6, 25, -11, 15]

  CH = cholfact(A)

  return CH \ B
end

function jacobi{T}(A::Matrix{T}, B::Vector{T}, n::Int, tol::Float64 = 1e-10)
  x = zeros(T, length(B))

  for i = 1:n
    x_new = zeros(length(B))
    for i = 1:size(A)[1]
      s1 = A[i, 1:i-1] * x[1:i-1]
      s2 = A[i, i+1:end] * x[i+1:end]
      x_new[i] = (B[i] - s1[1] - s2[1]) / A[i, i]
    end

    if isapprox(x, x_new, rtol = tol)
      break
    end

    x = x_new
  end

  return x
end


function main()
  stockreturns = [0.065, 0.0265, -0.0593, -0.001, 0.0346]
  mktreturns = [0.055, -0.09, -0.041, 0.045, 0.022]

  allvals = rand(9, 8)

  y_vals = allvals[:, 1]
  x_vals = allvals[:, 2:end]

  # return calc_beta(stockreturns, mktreturns)
  return ols_estimate(x_vals, y_vals)
end
