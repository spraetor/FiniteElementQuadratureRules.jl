struct Point <: AbstractDomain{0} end
struct Line <: AbstractDomain{1} end
struct Triangle <: AbstractDomain{2} end
struct Quadrilateral <: AbstractDomain{2} end
struct Tetrahedron <: AbstractDomain{3} end
struct Hexahedron <: AbstractDomain{3} end
struct Prism <: AbstractDomain{3} end
struct Pyramid <: AbstractDomain{3} end

const AllDomains = Union{Point,Line,Triangle,Quadrilateral,Tetrahedron,Hexahedron,Prism,Pyramid}

struct UnknownDomain{D} <: AbstractDomain{D} end

dimension(::Type{Domain}) where {D,Domain<:AbstractDomain{D}} = D

import Base: Symbol
Symbol(::Type{Point}) = :point
Symbol(::Type{Line}) = :line
Symbol(::Type{Triangle}) = :triangle
Symbol(::Type{Quadrilateral}) = :quadrilateral
Symbol(::Type{Tetrahedron}) = :tetrahedron
Symbol(::Type{Hexahedron}) = :hexahedron
Symbol(::Type{Prism}) = :prism
Symbol(::Type{Pyramid}) = :pyramid
Symbol(::UnknownDomain{D}) where D = :unknown


function domain(::Val{D}, region::String) where D
  if region == "simplex"
    if D == 1
      return Line
    elseif D == 2
      return Triangle
    elseif D == 3
      return Tetrahedron
    else
      return UnknownDomain{D}
    end
  elseif region == "cube"
    if D == 1
      return Line
    elseif D == 2
      return Quadrilateral
    elseif D == 3
      return Hexahedron
    else
      return UnknownDomain{D}
    end
  elseif region == "prism"
    @assert D == 3
    return Prism
  elseif region == "pyramid"
    @assert D == 3
    return Pyramid
  else
    return UnknownDomain{D}
  end
end
