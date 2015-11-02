using Distributions

function IR_deriv1(notional::Float64, K::Float64, alpha::Float64, sigma::Float64, r::Float64, dT::Float64, N::Int64, M::Int64)
  L = zeros(N+2, N+2) # forward rates
  D = zeros(N+2, N+2) # discount factors

  dW = zeros(N) # discount factors
  FV = zeros(N+1) # future value payment
  FVPrime = zeros(N+1) # numeraire-rebased FV payment
  # V = zeros(M) # simulation payoff

  W = Normal() # normal distribution for random num gen
  rands = zeros(N)

  df_prod = 1.0
  drift_sum = 0.0
  sumPV = 0.0
  PV = 0.0

  # initialize rates
  for i = 1:N+1
    L[i,1] = r
  end

  # main MC loop
  for nsim = 1:M
    payoff = 0.0
    # fill rands
    rand!(W, rands)
    for i = 1:N
      dW[i] = sqrt(dT) * rands[i]
    end

    # compute forward rates
    for n = 1:N, i = n + 1:N + 1
      drift_sum = 0.0
      for k = i + 1:N + 1
        drift_sum += (alpha * sigma * L[k,n]) / (1 + alpha * L[k,n])
      end
      L[i, n + 1] = L[i, n] * exp((-drift_sum * sigma - 0.5 * sigma ^ 2) * dT + sigma * dW[n])
    end

    # compute discount rates
    for n = 1:N + 1, i = n + 1:N + 2
      df_prod = 1.0
      for k = n:i
        df_prod *= 1 / (1 + alpha * L[k,n])
      end
      D[i,n] = df_prod
    end

    # compute everything else
    for i = 1:N + 1
      FV[i] = notional * alpha * (L[i,i] - K)
      FVPrime[i] = FV[i] * D[i + 1,i] / D[N + 2, i]
      payoff += FVPrime[i] * D[i + 1, 1]
    end

    sumPV += payoff
  end # end MC loop

  PV = sumPV / M

  return PV
end
