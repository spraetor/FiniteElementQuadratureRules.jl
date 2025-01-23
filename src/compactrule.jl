"""
  CompactQuadratureRule{T,Ω}
"""
struct CompactQuadratureRule{T<:Real, Ω<:AbstractDomain}
  degree::Int
  orbits::Vector{Int}
  values::Vector{T}
end

ctype(::CompactQuadratureRule{T,Ω}) where {T<:Real,Ω<:AbstractDomain} = T
domaintype(::CompactQuadratureRule{T,Ω}) where {T<:Real,Ω<:AbstractDomain} = Ω

# expand a compact rule into a quadrature rule
function expand(cqr::CompactQuadratureRule{T,Triangle}) where {T}
  nothing
end