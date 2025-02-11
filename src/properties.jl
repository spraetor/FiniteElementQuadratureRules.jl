function getProperties(domain::Ω, points::Vector{SVector{dim,T}}, weights::Vector{T}) where {Ω<:AbstractDomain,dim,T<:Real}
  properties = Symbol[]
  if all(weights .> T(0))
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