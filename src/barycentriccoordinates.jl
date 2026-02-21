"""
  barycentriccoordinates(domain::AbstractDomain, x::AbstractVector)

Transform reference element coordinates in a given domain into
barycentric coordinates. This is in particular useful for simplex
domains where these two coordinates differ. Also for prisms and pyramids
where some coordinates parts are associated to a conical extension.


Parameters:
 - `domain::AbstractDomain`: The domain of the original coordinates
 - `x::AbstractVector`: The coordinates in the domain

Result: Barycentric coordinates associated to the domain.
"""
function barycentriccoordinates(::AbstractDomain, x::AbstractVector)
  return x
end

# Specialization for simplex domains
function barycentriccoordinates(domain::AbstractSimplex, x::AbstractVector)
  return StaticVector{dimension(domain),eltype(x)}(1-sum(x), x...)
end