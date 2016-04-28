using StatsFuns
using Gadfly

function d1f(St::Float64, K::Float64, t::Float64, T::Float64, r::Float64, sigma::Float64)
  # Black scholes merton d1 function
  d1 = (log(St / K) + (r + 0.5 * sigma ^ 2) * (T - t)) / (sigma * sqrt(T - t))

  return d1
end

function BSM_call_value(St::Float64, K::Float64, t::Float64, T::Float64, r::Float64, sigma::Float64)
  # Calculates Black-Scholes Merton european call option value

  # Parameters
  # St : float - stock/index level at time t
  # K  : float - strike price
  # t  : float - valuation date
  # T  : float - date of maturity / time to maturity if t = 0; T > t
  # r  : float - constant, risk-less short rate
  # sigma : float - volatility

  # Returns
  # call-value : float - European call present at time t

  d1 = d1f(St, K, t, T, r, sigma)
  d2 = d1 - sigma * sqrt(T - t)

  call_value = St * normcdf(d1) - exp(-r * (T - t)) * K * normcdf(d2)

  return call_value
end

function main()
  K = 8000.0 # strike price
  T = 1.0 # Time to maturity
  r = 0.025 # constant risk-free short rate
  vol = 0.2 # constant volatility

  # Sample data generation
  S = linspace(4000, 12000, 150) # vector of index level values
  h = max(S - K, 0.0) # inner value
  C = [BSM_call_value(S0, K, 0.0, T, r, vol) for S0 in S]
  plot(x=S, y=C, Geom.line)
end
