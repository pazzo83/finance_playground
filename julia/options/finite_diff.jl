# for interpolation
using Dierckx

type FiniteDifference
  S0::Float64
  K::Float64
  r::Float64
  T::Float64
  sigma::Float64
  Smax::Float64
  M::Int64
  N::Int64
  is_call::Bool
  ds::Float64
  dt::Float64

  function FiniteDifference(S0::Float64, K::Float64, r::Float64, T::Float64, sigma::Float64, Smax::Float64, M::Int64, N::Int64; is_call::Bool = true)
    ds = Smax / M
    dt = T / N

    new(S0, K, r, T, sigma, Smax, M, N, is_call, ds, dt)
  end
end

# initial set up
function fd_initial_setup(M::Int64, N::Int64, Smax::Float64)
  i_values = collect(0:M - 1)
  j_values = collect(0:N - 1)
  boundary_conds = collect(linspace(0, Smax, M + 1))

  return i_values, j_values, boundary_conds
end

# setup the payoff grid
function fd_setup_grid(fd::FiniteDifference, j_values::Array{Int64, 1}, boundary_conds::Array{Float64, 1})
  grid = zeros(fd.M + 1, fd.N + 1)

  #setup boundary conditions
  if fd.is_call
    # right boundary
    for i = 1:fd.M + 1
      grid[i, end] = max(boundary_conds[i] - fd.K, 0)
    end

    # bottom boundary
    for j = 1:fd.N
      grid[end, j] = (fd.Smax - fd.K) * exp(-fd.r * fd.dt * (fd.N - j_values[j]))
    end
  else
    # right boundary
    for i = 1:fd.M + 1
      grid[i, end] = max(fd.K - boundary_conds[i], 0)
    end

    #top boundary
    for j = 1:fd.N
      grid[1, j] = (fd.K - fd.Smax) * exp(-fd.r * fd.dt * (fd.N - j_values[j]))
    end
  end

  return grid
end

# pricing methods
function price(fd::FiniteDifference)
  i_values, j_values, boundary_conds = fd_initial_setup(fd.M, fd.N, fd.Smax)
  grid = fd_setup_grid(fd, j_values, boundary_conds)

  # setup coefficients
  a = zeros(fd.M)
  b = zeros(fd.M)
  c = zeros(fd.M)

  # build a b c coefficients
  for i = 1:fd.M
    a[i] = 0.5 * (fd.r * fd.dt * i_values[i] - (fd.sigma ^ 2) * fd.dt * (i_values[i] ^ 2))
    b[i] = 1 + (fd.sigma ^ 2) * fd.dt * (i_values[i] ^ 2) + fd.r * fd.dt
    c[i] = -0.5 * (fd.r * fd.dt * i_values[i] + (fd.sigma ^ 2) * fd.dt * (i_values[i] ^ 2))
  end

  # build the diagonals
  # coeffs = diagm(a[3:fd.M], -1) + diagm(b[2:fd.M]) + diagm(c[c:fd.M - 1], 1)
  coeffs = Tridiagonal(a[3:fd.M], b[2:fd.M], c[2:fd.M - 1])

  # use LU factorization to solve linear systems of equations
  Alu = lufact(coeffs)
  aux = zeros(fd.M - 1)

  for j = fd.N:-1:1
    aux[1] = -a[2] * grid[1, j]
    grid[2:fd.M, j] = Alu \ (grid[2:fd.M, j + 1] + aux)
  end

  # interpolation
  return Spline1D(boundary_conds, grid[:, 1])(fd.S0)
end

function price_cn(fd::FiniteDifference)
  i_values, j_values, boundary_conds = fd_initial_setup(fd.M, fd.N, fd.Smax)
  grid = fd_setup_grid(fd, j_values, boundary_conds)

  # setup coefficients
  alpha = zeros(fd.M)
  beta = zeros(fd.M)
  gamma = zeros(fd.M)

  # build alpha beta gamma coefficients
  for i = 1:fd.M
    alpha[i] = 0.25 * fd.dt * ((fd.sigma ^ 2) * (i_values[i] ^ 2) - fd.r * i_values[i])
    beta[i] = -fd.dt * 0.5 * ((fd.sigma ^ 2) * (i_values[i] ^ 2) + fd.r)
    gamma[i] = 0.25 * fd.dt * ((fd.sigma ^ 2) * (i_values[i] ^ 2) + fd.r * i_values[i])
  end

  # build diagonal grids
  # M1 = -diagm(alpha[3:fd.M], -1) + diagm(1 - beta[2:fd.M]) - diagm(gamma[2:fd.M - 1], 1)
  # M2 = diagm(alpha[3:fd.M], -1) + diagm(1 + beta[2:fd.M]) + diagm(gamma[2:fd.M - 1], 1)
  M1 = Tridiagonal(-alpha[3:fd.M], 1 - beta[2:fd.M], -gamma[2:fd.M - 1])
  M2 = Tridiagonal(alpha[3:fd.M], 1 + beta[2:fd.M], gamma[2:fd.M - 1])

  # LU factorization to solve linear systems of equations
  Alu = lufact(M1)

  for j = fd.N:-1:1
    grid[2:fd.M, j] = Alu \ (M2 * grid[2:fd.M, j + 1])
  end

  # interpolation
  return Spline1D(boundary_conds, grid[:, 1])(fd.S0)
end
