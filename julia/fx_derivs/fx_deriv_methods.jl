function FX_deriv_finite_diff(S0::Float64, T::Float64, K::Float64, sigma::Float64, r::Float64, dx::Float64, dt::Float64, N::Int64, M::Int64)
  dtau = dt * (0.5 * sigma ^2)
  alpha = dtau / (dx ^ 2)
  k = r / (0.5 * sigma ^ 2)

  xmin = -1.0
  x = zeros(N)
  S = zeros(N)

  t = zeros(M)
  tau = zeros(M)

  u = zeros(N, M)
  v = zeros(N, M)

  # setup mesh - x and tau
  for i=1:N
    x[i] = xmin + (i - 1) * dx
    S[i] = K * exp(x[i])
  end

  for j=1:M
    t[j] = (j - 1) * dt
    tau[j] = (T - t[j]) / (0.5 * sigma ^ 2)
  end

  # setup initial conditions
  for i=1:N
    u[i, 1] = max(exp(0.5 * (k + 1) * x[i]) - exp(0.5 * (k - 1) * x[i]), 0.0)
  end

  # setup boundaries
  for j=2:M
    u[1, j] = 0.0
    u[N, j] = u[N, 1]
  end

  # compute forward differences
  for j = 1:M-1, i=2:N-1
    u[i, j + 1] = alpha * u[i + 1, j] + (1 - 2 * alpha) * u[i, j] + alpha * u[i - 1, j]
  end

  # transform solution
  for j=1:M, i=1:N
    v[i, j] = K ^ (0.5 * (1 + k)) * S[i] ^ (0.5 * (1 - k)) * exp(1.0 / 8.0 * (k + 1)^2 * sigma^2 * (T - t[i])) * u[i, j]
  end

  return dtau, alpha, k, x, S, t, tau, u, v
end

# FX_deriv_finite_diff(75.0, 0.5, 75.0, 0.3, 0.05, 0.5, 0.1, 5, 6)
