function barycentriccoordinates(domain::Ω, x::AbstractVector) where {Ω<:AbstractDomain}
  return x
end

function barycentriccoordinates(domain::Ω, x::AbstractVector) where {Ω<:AbstractSimplex}
  return StaticVector{dimension(domain),eltype(x)}(1-sum(x), x...)
end