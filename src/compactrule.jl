"""
  CompactQuadratureRule{T,Ω}
"""
struct CompactQuadratureRule{Ω<:AbstractDomain, T<:Real}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  positions::Vector{T}
end

ctype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω

# expand a compact rule into a quadrature rule
function expand(cqr::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  j = 1
  Point = typeof(sos[1].orbit()[1])
  points = Point[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = so.args
    for _ in 1:orbits
      push!(points, so.orbit(cqr.positions[j:j+n-1]...)...)
      j = j + n
    end
  end

  # maybe compute the weights here directly
  weights = getWeights(T,cqr.domain,cqr.degree,points,cqr.orbits)
  QuadratureRule(cqr.domain, cqr.degree, points, weights)
end

"""
  CompactQuadratureRuleWithWeights{T,Ω}

Compact rule that stores also weights.
"""
struct CompactQuadratureRuleWithWeights{Ω<:AbstractDomain, T<:Real}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  positions::Vector{T}
  weights::Vector{T}
end

ctype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω

# expand a compact rule with weights into a full quadrature rule
function expand(cqr::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  j = 1
  k = 1
  Point = typeof(sos[1].orbit()[1])
  points = Point[]
  weights = T[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = so.args
    for _ in 1:orbits
      push!(points, so.orbit(cqr.positions[j:j+n-1]...)...)
      append!(weights, fill(cqr.weights[k],so.size))
      j = j + n
      k = k + 1
    end
  end

  @assert length(points) == length(weights)
  QuadratureRule(cqr.domain, cqr.degree, points, weights)
end