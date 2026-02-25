using StaticArrays: SVector
using Printf: @sprintf

"""
    QuadratureRule{Ω,T,Point}

A `QuadratureRule` on a given domain Ω is a collection of points {xᵢ} in the associated
reference element and weights {wᵢ}, such that for a given polynomial p the quadrature
formula ∑ᵢ p(xᵢ) wᵢ = ∫p(x) dΩ is exact, up to polynomials of a certain degree.
"""
struct QuadratureRule{Ω<:AbstractDomain, T<:Real, Point<:AbstractVector{T}}
  ref::ReferenceElement{Ω,T,Point}
  degree::Int
  points::Vector{Point}
  weights::Vector{T}
  properties::Vector{Symbol}
end

"""
    QuadratureRule(domain::AbstractDomain, degree::Integer, points::Vector, weights::Vector)

Construct a new `QuadratureRule` and compute the properties of the rule automatically.
"""
function QuadratureRule(domain::Ω, degree::Integer, points::Vector{Point}, weights::Vector{T}) where {Ω<:AbstractDomain, Point<:AbstractVector, T<:Real}
  ref = ReferenceElement(domain)
  properties = getProperties(ref,points,weights)
  QuadratureRule(ReferenceElement{Ω,T,Point}(ref),degree,points,weights,properties)
end

"""
    QuadratureRule(domain::AbstractDomain, degree::Integer, points::Vector)

Construct a new `QuadratureRule` and compute the weights and properties of the
rule automatically.
"""
function QuadratureRule(domain::Ω, degree::Integer, points::Vector{Point}) where {Ω<:AbstractDomain, Point<:AbstractVector}
  ref = ReferenceElement(domain)
  weights = getWeights(ref,degree,points)
  properties = getProperties(ref,points,weights)
  QuadratureRule(ReferenceElement{Ω,eltype(weights),Point}(ref),degree,points,weights,properties)
end

# The type of the coordinates
ctype(::QuadratureRule{Ω,T,P}) where {Ω<:AbstractDomain,T,P} = T

# The dimension of the domain the quadrature rule is defined in
dimension(qr::QuadratureRule{Ω}) where {Ω<:AbstractDomain} = dimension(Ω)

# The type of the domain the quadrature rule is defined in
domaintype(::QuadratureRule{Ω}) where {Ω<:AbstractDomain} = Ω
domain(::QuadratureRule{Ω}) where {Ω<:AbstractDomain} = Ω()

isPositive(qr::QuadratureRule) = :positive in qr.properties
isInside(qr::QuadratureRule) = :inside in qr.properties
isPI(qr::QuadratureRule) = isPositive(qr) && isInside(qr)

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

  dom::AbstractDomain = domain(data["dim"],data["region"])
  degree::Integer = data["degree"]
  points = [ point(T,D,coords) for coords in data["coordinates"] ]
  weights = [ _parse(T,w) for w in data["weights"] ]
  properties = [ Symbol(p) for p in data["properties"] ]

  ref = ReferenceElement(dom)
  Point = eltype(points)
  QuadratureRule(ReferenceElement{domaintype(ref),T,Point}(ref),degree,points,weights,properties)
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
  println(io, "  domain = $(domain(qr)), ")
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
function Base.Dict(qr::QuadratureRule; reference::String="unknown", precision::Int=50)
  Dict(
    "reference" => reference,
    "region" => region(domain(qr)),
    "dim" => dimension(domain(qr)),
    "degree" => qr.degree,
    "properties" => String[ string(prop) for prop in qr.properties ],
    "coordinates" => [ String[ @sprintf("%0.*e", precision,pᵢ) for pᵢ in p ] for p in qr.points ],
    "weights" => String[ @sprintf("%0.*e", precision,w) for w in qr.weights ]
    )
end


function write_file(file::AbstractString, qr::QuadratureRule; reference::String="unknown", precision::Integer=50)
  open(file, "w") do f
    write(f, "reference: '$(reference)'\n")
    write(f, "region: $(region(domain(qr)))\n")
    write(f, "dim: $(dimension(domain(qr)))\n")
    write(f, "degree: $(qr.degree)\n")
    write(f, "properties: [$(length(qr.properties)>0 ? string(qr.properties[1]) : "")")
    for i in 2:length(qr.properties)
      write(f, ", $(string(qr.properties[i]))")
    end
    write(f, "]\n")
    write(f, "coordinates:\n")
    for p in qr.points
      write(f, "  - ['$(@sprintf("%0.*e",precision,p[1]))'")
      for i in 2:length(p)
        write(f, ", '$(@sprintf("%0.*e",precision,p[i]))'")
      end
      write(f, "]\n")
    end
    write(f, "weights:\n")
    for w in qr.weights
      write(f, "  - '$(@sprintf("%0.*e",precision,w))'\n")
    end
  end
end
