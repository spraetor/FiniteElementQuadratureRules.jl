module FiniteElementQuadratureRules

abstract type AbstractLocalBasis end
abstract type AbstractGeometry end
abstract type AbstractDomain end
abstract type AbstractSimplex <: AbstractDomain end
abstract type AbstractCube <: AbstractDomain end


import Base: parse
_parse(::Type{<:AbstractString}, x::AbstractString) = x
_parse(::Type{T}, x::AbstractString) where T<:Number = Base.parse(T,x)


include("domain.jl")
include("barycentriccoordinates.jl")
include("referenceelement.jl")
include("lagrange.jl")
include("geometry.jl")
include("symmetryorbits.jl")
include("polyset.jl")
include("pointset.jl")
include("weights.jl")
include("quadraturerule.jl")
include("compactrule.jl")
include("generate.jl")
include("integrate.jl")
include("transform.jl")
include("properties.jl")
include("dune.jl")

export AbstractCube, AbstractDomain, AbstractGeometry, AbstractLocalBasis,
  AbstractSimplex, AffineGeometry, BarycentricMonomial, CompactQuadratureRule, CompactQuadratureRuleWithWeights,
  Hexahedron, LagrangeLocalBasis, Line, MultiLinearGeometry, Point, PolySet, Prism, Pyramid,
  QuadratureRule, Quadrilateral, ReferenceElement, SymmetryOrbit, Tetrahedron, Triangle
export ctype, dimension, domain, domaintype, duneReferenceElement, expand, expandall, facets, generate,
  getProperties, getTuples, getWeights, integrate, order, position, symmetryOrbits, transform, vertices, volume,
  write_file

export LagrangeBasis

end # end module FiniteElementQuadratureRules
