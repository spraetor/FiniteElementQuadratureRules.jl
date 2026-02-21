struct Point <: AbstractSimplex end
struct Line <: AbstractCube end
struct Triangle <: AbstractSimplex end
struct Quadrilateral <: AbstractCube end
struct Tetrahedron <: AbstractSimplex end
struct Hexahedron <: AbstractCube end
struct Prism <: AbstractDomain end
struct Pyramid <: AbstractDomain end

const AllDomains = Union{Point,Line,Triangle,Quadrilateral,Tetrahedron,Hexahedron,Prism,Pyramid}

struct UnknownDomain <: AbstractDomain end

import Base: Symbol
Symbol(::Type{Point}) = :point
Symbol(::Type{Line}) = :line
Symbol(::Type{Triangle}) = :triangle
Symbol(::Type{Quadrilateral}) = :quadrilateral
Symbol(::Type{Tetrahedron}) = :tetrahedron
Symbol(::Type{Hexahedron}) = :hexahedron
Symbol(::Type{Prism}) = :prism
Symbol(::Type{Pyramid}) = :pyramid
Symbol(::Type{UnknownDomain}) = :unknown
Symbol(::Ω) where {Ω<:AbstractDomain} = Symbol(Ω)

dimension(::Type{Point}) = 0
dimension(::Type{Line}) = 1
dimension(::Type{Triangle}) = 2
dimension(::Type{Quadrilateral}) = 2
dimension(::Type{Tetrahedron}) = 3
dimension(::Type{Hexahedron}) = 3
dimension(::Type{Prism}) = 3
dimension(::Type{Pyramid}) = 3
dimension(::Type{UnknownDomain}) = 0
dimension(::Ω) where {Ω<:AbstractDomain} = dimension(Ω)

vertices(::Type{Point}) = 1
vertices(::Type{Line}) = 2
vertices(::Type{Triangle}) = 3
vertices(::Type{Quadrilateral}) = 4
vertices(::Type{Tetrahedron}) = 4
vertices(::Type{Hexahedron}) = 8
vertices(::Type{Prism}) = 6
vertices(::Type{Pyramid}) = 5
vertices(::Type{UnknownDomain}) = 0
vertices(::Ω) where {Ω<:AbstractDomain} = vertices(Ω)

facets(::Type{Point}) = 0
facets(::Type{Line}) = 2
facets(::Type{Triangle}) = 3
facets(::Type{Quadrilateral}) = 4
facets(::Type{Tetrahedron}) = 4
facets(::Type{Hexahedron}) = 6
facets(::Type{Prism}) = 5
facets(::Type{Pyramid}) = 5
facets(::Type{UnknownDomain}) = 0
facets(::Ω) where {Ω<:AbstractDomain} = facets(Ω)


# map a dim+region to a domain
function domaintype(dim::Int, region::Symbol)::Type{<:AbstractDomain}
  if dim == 0
    return Point
  elseif dim == 1
    return Line
  else
    if region == :simplex
      if dim == 2
        return Triangle
      elseif dim == 3
        return Tetrahedron
      else
        return UnknownDomain
      end
    elseif region == :cube
      if dim == 2
        return Quadrilateral
      elseif dim == 3
        return Hexahedron
      else
        return UnknownDomain
      end
    elseif region == :prism
      @assert dim == 3
      return Prism
    elseif region == :pyramid
      @assert dim == 3
      return Pyramid
    else
      return UnknownDomain
    end
  end
end

# map a dim+region to a domain, where the region might be given as a string
function domaintype(dim::Int, region::AbstractString)::Type{<:AbstractDomain}
  domaintype(dim, Symbol(region))
end

# construct a domain from dim + region
domain(dim, region)::AbstractDomain = domaintype(dim,region)()
