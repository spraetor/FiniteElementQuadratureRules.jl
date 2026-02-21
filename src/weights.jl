function getWeights(::Type{T}, domain::Ω, degree::Integer, points::AbstractVector{P}) where {T<:Real,Ω<:AbstractDomain,P<:AbstractVector}

  polyset = JacobiPolySet(domain, degree)
  A = zeros(T, length(polyset.basis), length(points))
  b = zeros(T, length(polyset.basis))

  for i in eachindex(polyset.basis)
    f = polyset.basis[i]
    for j in eachindex(points)
      A[i,j] = f(points[j])
    end
    b[i] = T(polyset.integrals[i])
  end

  weights = A\b
  return weights
end

function getWeights(domain::Ω, degree::Integer, points::AbstractVector{P}) where {Ω<:AbstractDomain,P<:AbstractVector}
  T = eltype(P)
  getWeights(T, domain, degree, points)
end



function getWeights(::Type{T}, domain::AbstractDomain, degree::Integer, points::AbstractVector{<:AbstractVector}, orbits::AbstractVector{<:Integer}) where {T<:Real}

  polyset = JacobiPolySet(domain, degree)
  so = symmetryOrbits(T,domain)
  nDifferentWeights = sum(orbits)
  A = zeros(T, length(polyset.basis), nDifferentWeights)
  b = zeros(T, length(polyset.basis))

  for i in eachindex(polyset.basis)
    f = polyset.basis[i]
    j = 1
    n = 1
    for k in eachindex(orbits)    # types of symmetry orbits
      for l in 1:orbits[k]        # number of orbits of this type
        fPointsSum = T(0)
        for m in 1:length(so[k])  # length of the orbit
          fPointsSum += f(points[n])
          n = n+1
        end
        A[i,j] = fPointsSum
        j = j+1
      end
    end
    b[i] = T(polyset.integrals[i])
  end

  w = A\b

  weights = T[]
  j = 1
  for k in eachindex(orbits)    # types of symmetry orbits
    for l in 1:orbits[k]        # number of orbits of this type
      append!(weights, fill(w[j],length(so[k])))
      j = j+1
    end
  end
  return weights
end

function getWeights(domain::AbstractDomain, degree::Integer, points::AbstractVector{P}, orbits::AbstractVector{<:Integer}) where {P<:AbstractVector}
  T = eltype(P)
  getWeights(T, domain, degree, points, orbits)
end