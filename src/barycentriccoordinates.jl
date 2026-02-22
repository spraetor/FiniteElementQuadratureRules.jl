using StaticArrays: SVector

"""
    barycentricCoordinates(domain::AbstractDomain, x::AbstractVector)

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
