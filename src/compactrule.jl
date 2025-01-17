struct CompactQuadratureRule{T<:Real, Ω<:AbstractDomain}
  degree::Int
  orbits::Vector{Int}
  values::Vector{T}
end

ctype(::CompactQuadratureRule{T,Ω}) where {T,Ω<:AbstractDomain} = T
domaintype(::CompactQuadratureRule{T,Ω}) where {T,Ω<:AbstractDomain} = Ω

function expand(cqr::CompactQuadratureRule{T,Triangle}) where {T}
  # expand a compact rule into a quadrature rule
end