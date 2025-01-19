using BibInternal: Entry
using BibParser: parse_file
using StaticArrays: SVector

struct QuadratureRule{T, D, Ω<:AbstractDomain}
  degree::Int
  points::Vector{SVector{D,T}}
  weights::Vector{T}
  properties::Vector{Symbol}
  accuracy::T
  bib::Entry
end

ctype(::QuadratureRule{T,D,Ω}) where {T,D,Ω<:AbstractDomain} = T
dimension(::QuadratureRule{T,D,Ω}) where {T,D,Ω<:AbstractDomain} = D
domaintype(::QuadratureRule{T,D,Ω}) where {T,D,Ω<:AbstractDomain} = Ω

import Base: parse
parse(::Type{<:AbstractString}, x::AbstractString) = x

function point(::Type{T}, ::Val{D}, coords::Vector{S}) where {T,D,S}
  @assert length(coords) == D
  SVector{D,T}((parse(T,c) for c in coords))
end

function (QuadratureRule{T})(dim::Int, region::String, data::Dict; bibFilename::String="") where T
  D::Val = Val(dim)

  degree = data["degree"]
  points = [ point(T,D,coords) for coords in data["coordinates"] ]
  weights = [ parse(T,w) for w in data["weights"] ]
  properties = [ Symbol(p) for p in data["properties"] ]
  accuracy = parse(T, data["accuracy"])

  if length(bibFilename) > 0
    bibFile = parse_file(bibFilename)
    bib = bibFile[data["reference"]]
  else
    bib = Entry(data["reference"], Dict("_type" => "unknown"))
  end

  QuadratureRule{T,dim,domain(D,region)}(degree,points,weights,properties,accuracy,bib)
end

import Base: length
length(qr::QuadratureRule) = length(qr.points)

import Base: show
function show(io::IO, qr::QuadratureRule)
  print(io, "{")
  print(io, "degree = $(qr.degree), ")
  # println(io, "points  = $(qr.points), ")
  # println(io, "weights = $(qr.weights), ")
  print(io, "props = $(qr.properties), ")
  # println(io, "accuracy= $(qr.accuracy), ")
  print(io, "bib = $(qr.bib.id) (year=$(qr.bib.date.year))")
  print(io, "}")
end
