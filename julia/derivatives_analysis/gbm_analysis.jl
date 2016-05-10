using DataFrames
using TimeSeries
using HypothesisTests

function dN(x::Float64, mu::Float64, sigma::Float64)
  z = (x - mu) / sigma
  pdf = exp(-0.5 * x ^ 2) / sqrt(2 * pi * sigma^2)

  return pdf
end

# Simulate a Number of Years of Daily Stock Quotes

function simulate_gbm()
  # Model params
  S0 = 100.0 # Initial index level
  T = 10.0 # time horizon
  r = 0.05 # risk-less short rate
  vol = 0.2 # instantaneous volatility

  # simulation params
  rng = MersenneTwister(250000)

  all_dates = Date(2004, 9, 30):Date(2014, 9, 30)

  gbm_dates = all_dates[Dates.dayofweek(all_dates) .< 6]

  M = length(gbm_dates) # Time steps
  I = 1 # index level paths
  dt = 1 / 252 # fixed for simplicity
  df = exp(-r * dt) # discount factor

  # stock price paths
  rands = randn(M, I)
  S = zeros(rands)
  S[1] = S0
  for t = 2:M
    S[t] = S[t-1] * exp((r - vol ^ 2 / 2) * dt + vol * rands[t] * sqrt(dt))
  end

  # gbm = DataFrame(Date=gbm_dates, index=S[:, 1])
  gbmTimeSeries = TimeArray(gbm_dates, S[:, 1], ["index"])

  rets = log(gbmTimeSeries["index"] ./ lag(gbmTimeSeries["index"], padding=true))

  gbm = DataFrame(Date=gbm_dates)
  gbm[:returns] = rets.values
  gbm[isnan(gbm[:returns]), :returns] = 0.0
  # Realized volatility (e.g. as defined for variance swaps)
  gbm[:rea_var] = 252 * cumsum(gbm[:returns] .^ 2) ./ range(1, M)
  gbm[:rea_vol] = sqrt(gbm[:rea_var])

  return gbm
end

function print_stats(data::DataFrame)
  ht = OneSampleTTest(data[:returns]) # for p value
  println("RETURN SAMPLE STATISTICS")
  println("----------------------------------------")
  println(@sprintf("Mean of Daily  Log Returns %9.6f", mean(data[:returns])))
  println(@sprintf("Std  of Daily  Log Returns %9.6f", std(data[:returns])))
  println(@sprintf("Mean of Annua. Log Returns %9.6f", mean(data[:returns]) * 252))
  println(@sprintf("Std  of Annua. Log Returns %9.6f", std(data[:returns]) * sqrt(252)))
  println("----------------------------------------")
  println(@sprintf("Skew of Sample Log Returns %9.6f", mean(data[:returns])))
  println(@sprintf("Kurt of Sample Log Returns %9.6f", kurtosis(data[:returns])))
  println("----------------------------------------")
  println(@sprintf("Normal test p-value        %9.6f", pvalue(ht)))
  println("----------------------------------------")
  println(@sprintf("Realized Volatility        %9.6f", data[:rea_vol][end]))
  println(@sprintf("Realized Variance          %9.6f", data[:rea_var][end]))
end
