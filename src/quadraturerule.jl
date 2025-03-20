using StaticArrays: SVector
import YAML: write_file

"""
  QuadratureRule{T,D,Ω}
"""
struct QuadratureRule{Ω<:AbstractDomain, T<:Real, Point<:AbstractVector{T}}
  domain::Ω
  degree::Int
  points::Vector{Point}
  weights::Vector{T}
  properties::Vector{Symbol}
end

# Construct a new QuadratureRule and compute the properties of the rule
function QuadratureRule(domain::Ω, degree::Integer, points::Vector{Point}, weights::Vector{T}) where {Ω<:AbstractDomain, Point<:AbstractVector, T<:Real}
  properties = getProperties(domain,points,weights)
  QuadratureRule(domain,degree,points,weights,properties)
end

# Construct a new QuadratureRule and compute the weights and properties of the rule
function QuadratureRule(domain::Ω, degree::Integer, points::Vector{Point}) where {Ω<:AbstractDomain, Point<:AbstractVector}
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

# convert a vector of coordinates represented as number or string into a static vector of type T
function point(::Type{T}, ::Val{D}, coords::Vector{S}) where {T,D,S}
  @assert length(coords) == D
  SVector{D,T}((_parse(T,c) for c in coords))
end

# Create a quadrature rules read from a json file
function (QuadratureRule{T})(dim::Integer, region::AbstractString, data::Dict) where T
  D::Val = Val(dim)

  domain = domain(D,region)
  degree = data["degree"]
  points = [ point(T,D,coords) for coords in data["coordinates"] ]
  weights = [ _parse(T,w) for w in data["weights"] ]
  properties = [ Symbol(p) for p in data["properties"] ]

  QuadratureRule(domain,degree,points,weights,properties)
end

# the length of a quadrature rule is equal to the number of quadrature points
import Base: length
length(qr::QuadratureRule) = length(qr.points)

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

function write_file(file::AbstractString, qr::QuadratureRule, data::Dict)
  qr_data = Dict(
    "reference" => data["reference"],
    "region" => data["region"],
    "dim" => data["dim"],
    "degree" => data["degree"],
    "properties" => String[ String(prop) for prop in qp.properties ],
    "coordinates" => [ String[ tostring(pᵢ) for pᵢ in p ] for p in qr.points ],
    "weights" => String[ tostring(w) for w in qr.weights ]
  )
  YAML.write_file(file, qr_data)
end