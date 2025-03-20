struct Point <: AbstractSimplex end
struct Line <: AbstractSimplex end
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

dimension(::Point) = 0
dimension(::Type{Point}) = 0
dimension(::Line) = 1
dimension(::Type{Line}) = 1
dimension(::Triangle) = 2
dimension(::Type{Triangle}) = 2
dimension(::Quadrilateral) = 2
dimension(::Type{Quadrilateral}) = 2
dimension(::Tetrahedron) = 3
dimension(::Type{Tetrahedron}) = 3
dimension(::Hexahedron) = 3
dimension(::Type{Hexahedron}) = 3
dimension(::Prism) = 3
dimension(::Type{Prism}) = 3
dimension(::Pyramid) = 3
dimension(::Type{Pyramid}) = 3
dimension(::UnknownDomain) = 0
dimension(::Type{UnknownDomain}) = 0

vertices(::Point) = 1
vertices(::Line) = 2
vertices(::Triangle) = 3
vertices(::Quadrilateral) = 4
vertices(::Tetrahedron) = 4
vertices(::Hexahedron) = 8
vertices(::Prism) = 6
vertices(::Pyramid) = 5
vertices(::UnknownDomain) = 0

facets(::Point) = 0
facets(::Line) = 2
facets(::Triangle) = 3
facets(::Quadrilateral) = 4
facets(::Tetrahedron) = 4
facets(::Hexahedron) = 6
facets(::Prism) = 5
facets(::Pyramid) = 5
facets(::UnknownDomain) = 0



function domaintype(dim::Int, region::Symbol)
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

function domaintype(dim::Int, region::AbstractString)
  domaintype(dim, Symbol(region))
end

domain(dim, region) = domaintype(dim,region)()