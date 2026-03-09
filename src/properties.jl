"""
    getProperties(domain::AbstractDomain, points::AbstractVector, weights::AbstractVector)

Collect properties of a quadrature given as vector of points and weights.
Currently, three properties are checked:
1. positive weights: property `:positive`
2. points strictly inside the domain: property `:inside`
3. points are inside or on the boundary of the domain: property `:boundary`

The function returns a vector of properties as corresponding Symbols.
"""
function getProperties(ref::ReferenceElement, points::AbstractVector, weights::AbstractVector)
  T = eltype(weights)
  properties = Symbol[]
  if all(weights .> zero(T))
    push!(properties, :positive)
  end

  if all((checkStrictlyInside(ref,p) for p in points))
    push!(properties, :inside)
  elseif all((checkInside(ref,p) for p in points))
    push!(properties, :boundary)
  end

  return properties
end


"""
    getProperties(qr::QuadratureRule)

Collect properties of a quadrature rule.
"""
getProperties(qr::QuadratureRule) = getProperties(qr.ref,qr.points,qr.weights)


"""
    getQuality(properties::AbstractVector{Symbol})
    getQuality(ref::ReferenceElement, points::AbstractVector, weights::AbstractVector)
    getQuality(qr::QuadratureRule)

Return the short quality marker associated with a quadrature rule.

For real-valued rules, the following markers are used:
- `:PI`: all weights are positive and all points are strictly inside the domain
- `:PB`: all weights are positive and all points are inside or on the boundary
- `:PO`: all weights are positive and at least one point is outside the domain
- `:NI`: at least one weight is non-positive and all points are strictly inside
- `:NB`: at least one weight is non-positive and all points are inside or on the boundary
- `:NO`: at least one weight is non-positive and at least one point is outside the domain

The literature also uses `:PC`, `:NC`, and `:CC` for rules with complex-valued
coordinates or weights. These markers are not produced by `QuadratureRule`,
which stores real-valued rules only.
"""
function getQuality(properties::AbstractVector{Symbol})
  prefix = :positive in properties ? "P" : "N"
  suffix = :inside in properties ? "I" : :boundary in properties ? "B" : "O"
  Symbol(prefix * suffix)
end

getQuality(ref::ReferenceElement, points::AbstractVector, weights::AbstractVector) =
  getQuality(getProperties(ref, points, weights))

getQuality(qr::QuadratureRule) = getQuality(qr.properties)
