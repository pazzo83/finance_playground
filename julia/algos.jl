# Algos

# Constants

const M_PERM = [3, 3, 1, 6, 4, 6, 8, 5]

function ins_sort_rec!{T <: Real}(seq::Vector{T}, i::Int)
  if i == 1
    return # base case, do nothing
  end

  ins_sort_rec!(seq, i-1) # sort 1...i

  j = i

  while j > 1 && seq[j-1] > seq[j]
    seq[j-1], seq[j] = seq[j], seq[j-1]
    j -= 1
  end

  return seq
end

ins_sort_rec{T <: Real}(seq::Vector{T}) = ins_sort_rec!(copy(seq), length(seq))

function ins_sort!{T <: Real}(seq::Vector{T})
  for i = 2:length(seq)
    j = i
    while j > 1 && seq[j-1] > seq[j]
      seq[j-1], seq[j] = seq[j], seq[j-1]
      j -= 1
    end
  end

  return seq
end

ins_sort{T <: Real}(seq::Vector{T}) = ins_sort!(copy(seq))

function sel_sort_rec!{T <: Real}(seq::Vector{T}, i::Int)
  if i == 1
    return # base case, do nothing
  end

  max_j = i # idx of largest value so far
  for j = 1:i # look for a larger value
    if seq[j] > seq[max_j]
      max_j = j # found one
    end
  end

  seq[i], seq[max_j] = seq[max_j], seq[i] # switch largest into place
  sel_sort_rec!(seq, i-1)

  return seq
end

sel_sort_rec{T <: Real}(seq::Vector{T}) = sel_sort_rec!(copy(seq), length(seq))

function sel_sort!{T <: Real}(seq::Vector{T})
  for i = length(seq):-1:1
    max_j = i
    for j = 1:i
      if seq[j] > seq[max_j]
        max_j = j
      end

      seq[i], seq[max_j] = seq[max_j], seq[i]
    end
  end

  return seq
end

sel_sort{T <: Real}(seq::Vector{T}) = sel_sort!(copy(seq))

function naive_max_perm(M::Vector{Int}, A::Set{Int} = Set{Int}())
  if isempty(A) # is the remaining people set not supplied
    A = Set(1:length(M)) # A = {1, 2, .., n}
  end

  if length(A) == 1 # base case - single person A
    return A
  end

  B = Set(Int[M[i] for i in A]) # The "pointed to" elements
  C = setdiff(A, B) # the "not pointed to" elements
  if ~isempty(C) # Any useless elements
    delete!(A, pop!(C)) # remove one
    return naive_max_perm(M, A) # solve remaining problem
  end

  return A # all usefull, return all
end

const G = rand(0:1, 100, 100)

function naive_celeb(G::Matrix{Int})
  n, k = size(G)

  for i = 1:k
    retVal = 0
    for j = 1:n
      if j == i
        continue
      end

      if G[i, j] > 0
        retVal = -1
        break
      end

      if G[j, i] == 0
        retVal = -1
        break
      end
    end

    if retVal == 0
      return i
    end
  end

  return -1
end

function celeb(G::Matrix{Int})
  n = size(G)[1]

  u, v = 1, 2 # The first two

  for c = 3:n+1 # Others to check
    if G[u, v] > 0 # u knows v?
      u = c # replace u
    else
      v = c # replace v
    end
  end

  if u == n+1 # u was replaced last
    c = v # use v
  else
    c = u # otherwise, u i s a candidate
  end

  retVal = 0
  for v = 1:n # for everyone else...
    if c == v # same person?
      continue # Skip
    end

    if G[c, v] > 0 # candidate knows other
      retVal = -1
      break
    end

    if G[v, c] == 0 # Other doesn't know candidate
      retVal = -1
      break
    end
  end

  if retVal == -1
    return -1
  else
    return c
  end
end

function naive_topsort(G::Dict{Int, IntSet}, S::IntSet = IntSet())
  if isempty(S)
    S = IntSet(keys(G)) # Default: ALL nodes
  end

  if length(S) == 1
    return collect(S) # Base case, single node
  end

  v = pop!(S) # Reduction: remove a node
  seq = naive_topsort(G, S) # recursion (assumption), n - 1
  min_i = 1
  for i in eachindex(seq)
    if v in G[seq[i]]
      min_i = i + 1 # after all dependencies
    end
  end

  insert!(seq, min_i, v)
  return seq
end

function topsort{T}(G::Dict{T, Set{T}})
  count = Dict{T, Int}([u => 0 for u in keys(G)]) # The in-degree for each node
  for u in keys(G)
    for v in G[u]
      count[v] += 1 # Count every valid in-edge
    end
  end

  count_zero(i::T, ::Set{T}) = count[i] == 0
  Q = T[u for (u, i) in filter(count_zero, G)] # Valid initial nodes
  S = T[] # result array

  while ~isempty(Q) # while we have start nodes...
    u = pop!(Q) # pick one
    push!(S, u) # use it as the first of the rest of our result
    for v in G[u]
      count[v] -= 1 # "uncount" its out-edges
      if count[v] == 0 # new valid start node?
        push!(Q, v) # Deal with them next
      end
    end
  end

  return S
end
