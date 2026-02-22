"""
    CompactQuadratureRule{T,Ω}

A compact (quadrature) rule is defined in terms of a vector of symmetry orbits
and their corresponding arguments, which allows to generate a list of quadrature
points in the reference element and also their corresponding quadrature weights.

The rule is defined on a domain, e.g., a triangle or quad, has an expected degree,
which means the maximal polynomial degree up to which the quadrature rule is exact.
"""
struct CompactQuadratureRule{Ω<:AbstractDomain, T<:Real}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  positions::Vector{T}
end


"""
    CompactQuadratureRule(::Type{T}, data::Dict)

Construct a `CompactQuadratureRuleWithWeights` from a YAML/Dict of strings and
string arrays. This is a convenience constructor typically used when reading a
quadrature rule from a YAML file. All information is encoded in the Dict, in
particular the fields
- `dim` and `region`: characterizing the domain
- `degree`: the quadrature degree
- `orbits`: the symmetry orbits
- `positions` or `arguments`: representing the arguments to the symmetry orbits
"""
function CompactQuadratureRule(::Type{T}, data::Dict) where T<:Real
  dom = domain(data["dim"],data["region"])
  degree = Int(data["degree"])
  orbits = Int[ o for o in data["orbits"] ]
  if haskey(data, "positions")
    positions = isnothing(data["positions"]) ? T[] : T[ _parse(T,p) for p in data["positions"] ]
  elseif haskey(data, "arguments")
    positions = isnothing(data["arguments"]) ? T[] : T[ _parse(T,p) for p in data["arguments"] ]
  else
    return nothing
  end

  CompactQuadratureRule{typeof(dom),T}(dom,degree,orbits,positions)
end


"""
    CompactQuadratureRule(data::Dict)

Construct a `CompactQuadratureRule` from parsed YAML/Dict data using `Float64`.
"""
CompactQuadratureRule(data::Dict) = CompactQuadratureRule(Float64, data)

ctype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω


"""
    expand(cqr::CompactQuadratureRule{Ω,T})

Expand a compact rule into a full quadrature rule. This first expands the symmetry
orbits to generate a list of quadrature points and then computes the associated
quadrature weights.
"""
function expand(cqr::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  if length(cqr.orbits) > length(sos)
    error("Number of orbits incompatible with available symmetric orbits.")
    return nothing
  end

  j = 1
  PointRef = SVector{dimension(Ω),T}
  Point = typeof(barycentricCoordinates(Ω(), zero(PointRef)))
  points = Point[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = args(so)
    for _ in 1:orbits
      append!(points, expand(so,view(cqr.positions,j:j+n-1)))
      j = j + n
    end
  end

  # maybe compute the weights here directly
  coords = transformCoordinates(cqr.domain,points)
  weights = getWeights(T,cqr.domain,cqr.degree,coords,cqr.orbits)
  QuadratureRule(cqr.domain, cqr.degree, coords, weights)
end


"""
    CompactQuadratureRuleWithWeights{T,Ω}

A compact (quadrature) rule with weights is defined in terms of a vector of
symmetry orbits, their corresponding arguments, and a list of associated weights.
This allows to generate a list of quadrature points in the reference element.

The rule is defined on a domain, e.g., a triangle or quad, has an expected degree,
which means the maximal polynomial degree up to which the quadrature rule is exact.
"""
struct CompactQuadratureRuleWithWeights{Ω<:AbstractDomain, T<:Real}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  positions::Vector{T}
  weights::Vector{T}
end


"""
  CompactQuadratureRuleWithWeights(::Type{T}, data::Dict)

Construct a `CompactQuadratureRuleWithWeights` from a YAML/Dict of strings and
string arrays. This is a convenience constructor typically used when reading a
quadrature rule from a YAML file. All information is encoded in the Dict, in
particular the fields
- `dim` and `region`: characterizing the domain
- `degree`: the quadrature degree
- `orbits`: the symmetry orbits
- `positions`: representing the arguments to the symmetry orbits
- `weights`: for the quadrature weights
"""
function CompactQuadratureRuleWithWeights(::Type{T}, data::Dict) where T
  dom = domain(data["dim"],data["region"])
  degree = Int(data["degree"])
  orbits = Int[ o for o in data["orbits"] ]
  positions = isnothing(data["positions"]) ? T[] : T[ _parse(T,p) for p in data["positions"] ]
  weights = isnothing(data["weights"]) ? T[] : T[ _parse(T,w) for w in data["weights"] ]

  CompactQuadratureRuleWithWeights(dom,degree,orbits,positions,weights)
end


"""
    CompactQuadratureRuleWithWeights(data::Dict)

Construct a `CompactQuadratureRuleWithWeights` from parsed YAML/Dict data using `Float64`.
"""
CompactQuadratureRuleWithWeights(data::Dict) = CompactQuadratureRuleWithWeights(Float64,data)

ctype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω


"""
    expand(cqr::CompactQuadratureRuleWithWeights{Ω,T})

Expand a compact rule into a full quadrature rule. This first expands the symmetry
orbits to generate a list of quadrature points and combines it with the given
quadrature weights.
"""
function expand(cqr::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  if length(cqr.orbits) > length(sos)
    return nothing
  end
  j = 1
  k = 1
  PointRef = SVector{dimension(Ω),T}
  Point = typeof(barycentricCoordinates(Ω(), zero(PointRef)))
  points = Point[]
  weights = T[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = args(so)
    for _ in 1:orbits
      append!(points, expand(so,view(cqr.positions,j:j+n-1)))
      append!(weights, fill(cqr.weights[k],length(so)))
      j = j + n
      k = k + 1
    end
  end

  @assert length(points) == length(weights)
  QuadratureRule(cqr.domain, cqr.degree,
    transformCoordinates(cqr.domain,points), transformWeights(cqr.domain, weights))
end