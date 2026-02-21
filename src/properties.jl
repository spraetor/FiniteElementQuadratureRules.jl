"""
  getProperties(domain::AbstractDomain, points::AbstractVector, weights::AbstractVector)

Collect properties of a quadrature given as vector of points and weights.
Currently, three properties are checked:
1. positive weights: property `:positive`
2. points strictly inside the domain: property `:inside`
3. points are inside or on the boundary of the domain: property `:boundary`

The function returns a vector of properties as corresponding Symbols.
"""
function getProperties(domain::AbstractDomain, points::AbstractVector, weights::AbstractVector)
  T = eltype(weights)
  properties = Symbol[]
  if all(weights .> zero(T))
    push!(properties, :positive)
  end

  ref = ReferenceElement(domain)
  if all((checkStrictlyInside(ref,p) for p in points))
    push!(properties, :inside)
  elseif all((checkInside(ref,p) for p in points))
    push!(properties, :boundary)
  end

  return properties
end