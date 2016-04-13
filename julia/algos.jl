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
