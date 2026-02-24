"""
    getWeights(::Type{T}, ref::ReferenceElement, degree::Integer, points::AbstractVector)

Compute the quadrature weight associated to given quadrature `points` for quadrature
in a given reference element `ref`, such that polynomial up to `degree` are integrated
exactly. The parameter `T` represents the number type to be used for the computation.
"""
function getWeights(::Type{T}, ref::ReferenceElement, degree::Integer, points::AbstractVector) where {T<:Real}

  polyset = JacobiPolySet(domain(ref), degree)
  A = Matrix{T}(undef, length(polyset.basis), length(points))
  b = Vector{T}(undef, length(polyset.basis))

  for j in eachindex(points)
    for i in eachindex(polyset.basis)
      f = polyset.basis[i]
      A[i,j] = f(points[j])
    end
  end
  for i in eachindex(polyset.basis)
    b[i] = T(polyset.integrals[i])
  end

  weights = A\b
  return weights
  # volIn = volume(ReferenceElement(domain(ref)))
  # volOut = volume(ref)
  # return map(w -> w*volOut/volIn, weights)
end

function getWeights(::Type{T}, domain::AbstractDomain, degree::Integer, points::AbstractVector) where {T<:Real}
  getWeights(T, ReferenceElement(domain), degree, points)
end


"""
    getWeights(::Type{T}, domain::AbstractDomain, degree::Integer, points::AbstractVector)

Compute the quadrature weight using Float64 precision.
"""
function getWeights(domain::AbstractDomain, degree::Integer, points::AbstractVector{P}) where {P<:AbstractVector}
  T = eltype(P)
  getWeights(T, ReferenceElement(domain), degree, points)
end

function getWeights(ref::ReferenceElement, degree::Integer, points::AbstractVector)
  T = eltype(P)
  getWeights(T, ref, degree, points)
end

"""
    getWeights(::Type{T}, ref::ReferenceElement, degree::Integer, points::AbstractVector, orbits::AbstractVector)

Compute the quadrature weight associated to given quadrature `points` for quadrature
in a given reference element `ref`, such that polynomial up to `degree` are integrated
exactly. This overload takes into account that the points are associated to symmetry
orbits and thus points in the same orbit share a quadrature weight. This makes the
computation more stable and faster. The parameter `T` represents the number type to be
used for the computation.
"""
function getWeights(::Type{T}, ref::ReferenceElement, degree::Integer, points::AbstractVector, orbits::AbstractVector) where {T<:Real}

  polyset = JacobiPolySet(domain(ref), degree)
  so = symmetryOrbits(T,domain(ref))
  nDifferentWeights = sum(orbits)
  A = Matrix{T}(undef, length(polyset.basis), nDifferentWeights)
  b = Vector{T}(undef, length(polyset.basis))

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
  # volIn = volume(ReferenceElement(domain(ref)))
  # volOut = volume(ref)
  # return map(w -> w*volOut/volIn, weights)
end

function getWeights(::Type{T}, domain::AbstractDomain, degree::Integer, points::AbstractVector, orbits::AbstractVector) where {T<:Real}
  getWeights(T, ReferenceElement(domain), degree, points, orbits)
end

"""
    getWeights(::Type{T}, domain::AbstractDomain, degree::Integer, points::AbstractVector, orbits::AbstractVector)

Compute the quadrature weight using Float64 precision.
"""
function getWeights(domain::AbstractDomain, degree::Integer, points::AbstractVector{P}, orbits::AbstractVector) where {P<:AbstractVector}
  T = eltype(P)
  getWeights(T, ReferenceElement(domain), degree, points, orbits)
end

function getWeights(ref::ReferenceElement, degree::Integer, points::AbstractVector, orbits::AbstractVector)
  T = eltype(P)
  getWeights(T, ref, degree, points, orbits)
end