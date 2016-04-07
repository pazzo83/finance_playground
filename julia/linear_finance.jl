# for CAPM
calc_beta(stockreturns::Vector{Float64}, mktreturns::Vector{Float64}) = linreg(stockreturns, mktreturns)

# for Arbitrage Pricing Theory model
ols_estimate(x_vals::Matrix{Float64}, y_vals::Vector{Float64}) = linreg(x_vals, y_vals)


function main()
  stockreturns = [0.065, 0.0265, -0.0593, -0.001, 0.0346]
  mktreturns = [0.055, -0.09, -0.041, 0.045, 0.022]

  allvals = rand(9, 8)

  y_vals = allvals[:, 1]
  x_vals = allvals[:, 2:end]

  # return calc_beta(stockreturns, mktreturns)
  return ols_estimate(x_vals, y_vals)
end
