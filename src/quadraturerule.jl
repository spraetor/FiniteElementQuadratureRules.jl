using StaticArrays: SVector

"""
  QuadratureRule{Ω,T,Point}

A `QuadratureRule` on a given domain Ω is a collection of points {xᵢ} in the associated
reference element and weights {wᵢ}, such that for a given polynomial p the quadrature
formula ∑ᵢ p(xᵢ) wᵢ = ∫p(x) dΩ is exact, up to polynomials of a certain degree.
"""
struct QuadratureRule{Ω<:AbstractDomain, T<:Real, Point<:AbstractVector{T}}
  domain::Ω
  degree::Int
  points::Vector{Point}
  weights::Vector{T}
  properties::Vector{Symbol}
end

"""
  QuadratureRule(domain::AbstractDomain, degree::Integer, points::Vector, weights::Vector)

Construct a new `QuadratureRule` and compute the properties of the rule automatically.
"""
function QuadratureRule(domain::AbstractDomain, degree::Integer, points::Vector{Point}, weights::Vector{T}) where {Point<:AbstractVector, T<:Real}
  properties = getProperties(domain,points,weights)
  QuadratureRule(domain,degree,points,weights,properties)
end

"""
  QuadratureRule(domain::AbstractDomain, degree::Integer, points::Vector)

Construct a new `QuadratureRule` and compute the weights and properties of the
rule automatically.
"""
function QuadratureRule(domain::AbstractDomain, degree::Integer, points::Vector{Point}) where {Point<:AbstractVector}
  weights = getWeights(domain,degree,points)
  properties = getProperties(domain,points,weights)
  QuadratureRule(domain,degree,points,weights,properties)
end

# The type of the coordinates
ctype(::QuadratureRule{Ω,T,P}) where {Ω<:AbstractDomain,T,P} = T

# The dimension of the domain the quadrature rule is defined in
dimension(qr::QuadratureRule{Ω,T,P}) where {Ω<:AbstractDomain,T,P} = dimension(qr.domain)

# The type of the domain the quadrature rule is defined in
domaintype(::QuadratureRule{Ω,T,P}) where {Ω<:AbstractDomain,T,P} = Ω
domain(::QuadratureRule{Ω,T,P}) where {Ω<:AbstractDomain,T,P} = domain(Ω)

# convert a vector of coordinates represented as number or string into a static vector of type T
function point(::Type{T}, ::Val{D}, coords::Vector{S}) where {T,D,S}
  @assert length(coords) == D
  SVector{D,T}((_parse(T,c) for c in coords))
end


"""
  QuadratureRule(::Type{T}, data::Dict)

Construct a `QuadratureRule` from a YAML/Dict of strings and
string arrays. This is a convenience constructor typically used when reading a
quadrature rule from a YAML file. All information is encoded in the Dict, in
particular the fields
- `dim` and `region`: characterizing the domain
- `degree`: the quadrature degree
- `coordinates`: representing the quadrature points
- `weights`: for the quadrature weights
- `properties`: characterizing properties of the quadrature rule
"""
function QuadratureRule(::Type{T}, data::Dict) where T
  D::Val = Val(data["dim"])

  domain::AbstractDomain = domain(D,data["region"])
  degree::Integer = data["degree"]
  points = [ point(T,D,coords) for coords in data["coordinates"] ]
  weights = [ _parse(T,w) for w in data["weights"] ]
  properties = [ Symbol(p) for p in data["properties"] ]

  QuadratureRule(domain,degree,points,weights,properties)
end

"""
  QuadratureRule(data::Dict)

Construct a `QuadratureRule` from parsed YAML/Dict data using `Float64`.
"""
QuadratureRule(data::Dict) = QuadratureRule(Float64, data)

# the length of a quadrature rule is equal to the number of quadrature points
import Base: length
length(qr::QuadratureRule) = length(qr.points)

# Print the quadrature rule
import Base: show
function show(io::IO, qr::QuadratureRule)
  println(io, "{")
  println(io, "  domain = $(qr.domain), ")
  println(io, "  degree = $(qr.degree), ")
  println(io, "  length = $(length(qr)), ")
  println(io, "  points  = $(qr.points), ")
  println(io, "  weights = $(qr.weights), ")
  println(io, "  props = $(qr.properties), ")
  println(io, "}")
end

import Base: Dict

"""
  Dict(qr::QuadratureRule, ref::String = "unknown")

Convert the given `QuadratureRule` into a Dict for exporting into a YAML file.
The optional parameter `ref` refers to a bibtex key used to reference a publication
where the quadrature rule is extracted from.
"""
function Base.Dict(qr::QuadratureRule, ref::String = "unknown")
  Dict(
    "reference" => ref,
    "region" => region(qr.domain),
    "dim" => dimension(qr.domain),
    "degree" => qr.degree,
    "properties" => String[ string(prop) for prop in qr.properties ],
    "coordinates" => [ String[ string(pᵢ) for pᵢ in p ] for p in qr.points ],
    "weights" => String[ string(w) for w in qr.weights ]
    )
end
