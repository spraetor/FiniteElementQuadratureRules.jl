using StaticArrays: @SMatrix, @SVector, SMatrix, SVector

"""
    barycentricCoordinates(domain::AbstractDomain, x::AbstractVector)
    barycentricCoordinates(ref::ReferenceElement, x::AbstractVector)

Transform reference element coordinates in a given domain into
barycentric coordinates. This is in particular useful for simplex
domains where these two coordinates differ. Also for prisms and pyramids
where some coordinates parts are associated to a conical extension.

Parameters:
 - `domain::AbstractDomain`: The domain of the original coordinates
 - `x::AbstractVector`: The coordinates in the reference domain

Result: Barycentric coordinates associated to the domain.
"""
function barycentricCoordinates(::AbstractDomain, x::AbstractVector)
  return x
end

function barycentricCoordinates(ref::ReferenceElement, x::AbstractVector)
  return barycentricCoordinates(domain(ref),x)
end

function barycentricCoordinates(ref::ReferenceElement{Ω}, x::AbstractVector) where {Ω<:AbstractSimplex}
  dim = dimension(Ω)
  n = dim + 1
  @assert length(x) == dim
  T = float(promote_type(eltype(x), ctype(ref)))
  coords = coordinates(ref)
  A = SMatrix{n,n,T}(ntuple(n*n) do k
    i = (k - 1) % n + 1
    j = (k - 1) ÷ n + 1
    i <= dim ? T(coords[j][i]) : one(T)
  end)
  b = SVector{n,T}(ntuple(n) do i
    i <= dim ? T(x[i]) : one(T)
  end)
  return A \ b
end

# Specializations for triangle reference domains
function barycentricCoordinates(::Triangle, x::AbstractVector)
  @assert length(x) == 2
  T = eltype(x)
  return SVector{3,T}(-(x[1]+x[2])/2, (one(T)+x[1])/2, (one(T)+x[2])/2)
end

# Specializations for tetrahedron reference domains
function barycentricCoordinates(::Tetrahedron, x::AbstractVector)
  @assert length(x) == 3
  T = eltype(x)
  return SVector{4,T}(
    -(one(T)+x[1]+x[2]+x[3])/2,
    (one(T)+x[1])/2,
    (one(T)+x[2])/2,
    (one(T)+x[3])/2,
  )
end

# Specializations for prism reference domains
function barycentricCoordinates(::Prism, x::AbstractVector)
  @assert length(x) == 3
  T = eltype(x)
  return SVector{4,T}(
    -(x[1]+x[2])/2,
    (one(T)+x[1])/2,
    (one(T)+x[2])/2,
    x[3],
  )
end

function barycentricCoordinates(ref::ReferenceElement{Prism}, x::AbstractVector)
  @assert length(x) == 3
  T = float(promote_type(eltype(x), ctype(ref)))
  coords = _transform(SVector{3,T}, coordinates(ref))
  A = SMatrix{3,3,T}(
    coords[1][1], coords[1][2], one(T),
    coords[2][1], coords[2][2], one(T),
    coords[3][1], coords[3][2], one(T),
  )
  λ = A \ SVector{3,T}(x[1], x[2], one(T))
  z_bottom = (coords[1][3] + coords[2][3] + coords[3][3]) / T(3)
  z_top = (coords[4][3] + coords[5][3] + coords[6][3]) / T(3)
  z = (T(2) * T(x[3]) - (z_top + z_bottom)) / (z_top - z_bottom)
  return SVector{4,T}(λ[1], λ[2], λ[3], z)
end


# transform the coordinates from the internal (barycentric) to the reference element domain
transformCoordinates(::AbstractCube, X::AbstractVector{P}) where {P<:AbstractVector} = X

function transformCoordinates(::Triangle, X::AbstractVector{P}) where {P<:AbstractVector}
  let A = @SMatrix [-1 1 -1; -1 -1 1]
    map(λ -> A*λ, X)
  end
end

function transformCoordinates(::Tetrahedron, X::AbstractVector{P}) where {P<:AbstractVector}
  let A = @SMatrix [-1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1]
    map(λ -> A*λ, X)
  end
end

function transformCoordinates(::Prism, X::AbstractVector{P}) where {P<:AbstractVector}
  let A = SMatrix{3,4}(-1,-1,0, 1,-1,0, -1,1,0, 0,0,1),
      b = SVector{3}(0, 0, 0)
    map(x -> A*x + b, X)
  end
end
transformCoordinates(::Pyramid, X::AbstractVector{P}) where {P<:AbstractVector} = X

transformCoordinates(domain::AbstractDomain, X::NTuple{N}) where {N} = transformCoordinates(domain, SVector(X))


# transform the quadrature weights when changing the coordinates from barycentric to reference domain
transformWeights(::AbstractCube, W::AbstractVector{<:Real}) = W
transformWeights(::Triangle, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Tetrahedron, W::AbstractVector{<:Real}) = map(w -> 4*w/3, W)
transformWeights(::Prism, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Pyramid, W::AbstractVector{<:Real}) = W
