using Distributions

function EQ1(S0::Float64, T::Float64, K::Float64, sigma::Float64, r::Float64, N::Int64, M::Int64)
  sum_payoff = 0.0
  dt = T / N

  S = zeros(N + 1)
  rands = zeros(N)
  W = Normal()

  for j = 1:M
    S[1] = S0
    rand!(W, rands)
    for i=1:N
      S[i + 1] = S[i] * (1 + r * dt + sigma * sqrt(dt) * rands[i])
    end

    sum_payoff += max(S[N] - K, 0.0)
  end

  premium = exp(-r * T) * (sum_payoff / M)

  premium
end

function EQ2(S01::Float64, S02::Float64, T::Float64, sigma1::Float64, sigma2::Float64, r::Float64, rho::Float64, N::Int64, M::Int64)
  sum_payoff = 0.0
  dt = T / N

  S1 = zeros(N + 1) # equity 1
  S2 = zeros(N + 1) # equity 2
  rands1 = zeros(N) # rands for eq 1
  rands2 = zeros(N) # rands for eq 2
  W = Normal() # normal distribution to draw rands from

  for j = 1:M
    S1[1] = S01
    S2[1] = S02
    rand!(W, rands1)
    rand!(W, rands2)

    for i=1:N
      S1[i + 1] = S1[i] * (1 + r * dt + sigma1 * sqrt(dt) * rands1[i])
      S2[i + 1] = S2[i] * (1 + r * dt + sigma2 * sqrt(dt) * (rands1[i] * rho  + sqrt(1 - rho * rho) * rands2[i]))
    end

    sum_payoff += max(S1[N], S2[N])
  end

  premium = exp(-r * T) * (sum_payoff / M)

  premium
end
