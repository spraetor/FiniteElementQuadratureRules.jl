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
include("referenceelement.jl")
include("lagrange.jl")
include("geometry.jl")
include("symmetryorbits.jl")
include("polyset.jl")
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
export ctype, dimension, domain, domaintype, duneReferenceElement, expand, facets, generate,
  getProperties, getTuples, getWeights, integrate, order, position, symmetryOrbits, transform, vertices, volume

end # end module FiniteElementQuadratureRules
