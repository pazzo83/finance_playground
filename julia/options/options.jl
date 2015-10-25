type StockOption
  S0::Float64
  K::Float64
  r::Float64
  T::Float64
  N::Int64
  pu::Float64
  pd::Float64
  div::Float64
  sigma::Float64
  df::Float64
  dt::Float64
  is_call::Bool
  is_euro::Bool

  function StockOption(S0::Float64, K::Float64, r::Float64, T::Float64, N::Int64; pu = 0.0, pd = 0.0, div = 0.0, sigma = 0.0, is_call = true, is_euro = true)
    # calculate dt
    dt = T / N # single time step in years

    # calculate DF
    df = exp(-(r - div) * dt) # e^-rt

    new(S0, K, r, T, N, pu, pd, div, sigma, df, dt, is_call, is_euro)
  end
end

# define pricing methods
abstract BinomialMethod

type BinomialNormal <: BinomialMethod end
type BinomialCRR <: BinomialMethod end
type BinomialLR <: BinomialMethod end

# build binomial tree
# function build_tree(s_opt::StockOption)
#   u = 1 + s_opt.pu
#   d = 1 - s_opt.pd
#   qu = (exp((s_opt.r - s_opt.div) * s_opt.dt) - d) / (u - d)
#   qd = 1 - qu
#
#   # init array
#   stockTree = Vector{Float64}[ zeros(m) for m = 1:s_opt.N + 1]
#
#   for i = 1:s_opt.N + 1
#     for j = 1:i
#       stockTree[i][j] = s_opt.S0 * u ^ (j-1) * d ^ (i - j)
#     end
#   end
#
#   return stockTree
# end

# function build_tree_matrix(s_opt::StockOption)
#   u = 1 + s_opt.pu
#   d = 1 - s_opt.pd
#   qu = (exp((s_opt.r - s_opt.div) * s_opt.dt) - d) / (u - d)
#   qd = 1 - qu
#
#   # init Matrix
#   stockMatrix = zeros(s_opt.N + 1, s_opt.N + 1)
#
#   for i =1:size(stockMatrix)[1]
#     for j = 1:i
#       stockMatrix[j, i] = s_opt.S0 * u ^ (j-1) * d ^ (i - j)
#     end
#   end
#
#   return stockMatrix
# end

function check_early_exercise(payoff::Float64, ST::Float64, K::Float64, is_call::Bool)
  return max(payoff, is_call ? ST - K : K - ST)
end

function build_tree_matrix(S0::Float64, n::Int64, u::Float64, d::Float64)
  stockTree = zeros(n + 1, n + 1)

  for i = 1:n + 1, j = 1:n + 1
    stockTree[i, j] = S0 * u ^ (j-1) * d ^ (i - 1)
  end

  return stockTree
end

function option_value_tree_matrix(K::Float64, n::Int64, binom::Array, qu::Float64, qd::Float64, df::Float64, is_call::Bool, is_euro::Bool)
  payoffs = zeros(n + 1, n + 1)

  for i = 1:n + 1
    payoffs[n+1 - (i - 1), i] = max(0, is_call ? binom[n+1 - (i - 1), i] - K : K - binom[n+1 - (i - 1), i])
  end

  for i = n:-1:1, j = n - (i - 1):-1:1
    payoffs[i,j] = (qu * payoffs[i, j+1] + qd * payoffs[i+1, j]) * df

    if ! is_euro
      payoffs[i,j] = check_early_exercise(payoffs[i,j], binom[i,j], K, is_call)
    end
  end

  return payoffs
end

# diff setup params for diff methods
function setup_params(s::StockOption, m::BinomialMethod)
  u = 1 + s.pu
  d = 1 - s.pd
  qu = (exp((s.r - s.div) * s.dt) - d) / (u - d)
  qd = 1 - qu
  return u, d, qu, qd
end

function setup_params(s::StockOption, m::BinomialCRR)
  u = exp(s.sigma * sqrt(s.dt))
  d = 1.0 / u
  qu = (exp((s.r - s.div) * s.dt) - d) / (u - d)
  qd = 1.0 - qu

  return u, d, qu, qd
end

function setup_params(s::StockOption, m::BinomialLR)
  odd_N = isodd(s.N) ? s.N : s.N + 1
  d1 = (log(s.S0 / s.K) + ((s.r - s.div) + (s.sigma ^ 2) / 2.0) * s.T) / (s.sigma * sqrt(s.T))
  d2 = (log(s.S0 / s.K) + ((s.r - s.div) - (s.sigma ^ 2) / 2.0) * s.T) / (s.sigma * sqrt(s.T))

  pp_2_inversion(z::Float64, n::Int64) = 0.5 + copysign(1, z) * sqrt(0.25 - 0.25 * exp( -((z / (n + 1.0 / 3.0 + 0.1 / (n + 1.0))) ^ 2) * (n + 1.0 / 6.0)))
  pbar = pp_2_inversion(d1, odd_N)
  p = pp_2_inversion(d2, odd_N)
  u = 1.0 / s.df * pbar / p
  d = (1 / s.df - p * u) / (1 - p)
  qu = p
  qd = 1 - p

  return u, d, qu, qd
end

function price(s_opt::StockOption; method::BinomialMethod = BinomialNormal())
  u, d, qu, qd = setup_params(s_opt, method)

  binom_tree = build_tree_matrix(s_opt.S0, s_opt.N, u, d)
  payoffs = option_value_tree_matrix(s_opt.K, s_opt.N, binom_tree, qu, qd, s_opt.df, s_opt.is_call, s_opt.is_euro)

  return payoffs[1,1]
end
