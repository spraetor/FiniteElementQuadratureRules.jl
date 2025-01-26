"""
  CompactQuadratureRule{T,Ω}
"""
struct CompactQuadratureRule{T<:Real, Ω<:AbstractDomain}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  values::Vector{T}
end

ctype(::CompactQuadratureRule{T,Ω}) where {T<:Real,Ω<:AbstractDomain} = T
domaintype(::CompactQuadratureRule{T,Ω}) where {T<:Real,Ω<:AbstractDomain} = Ω

# expand a compact rule into a quadrature rule
function expand(cqr::CompactQuadratureRule{T,Ω}) where {T<:Real,Ω<:AbstractDomain}
  sos = symmetryOrbits(T,cqr.domain)
  pos = 1
  points = SVector{3,T}[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    for _ in 1:orbits
      append!(points, so.orbit(cqr.values[pos:pos+so.args-1]...)...)
      pos = pos+so.args
    end
  end
  # TODO: Compute weights
end