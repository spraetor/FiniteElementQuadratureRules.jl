module FiniteElementQuadratureRules

abstract type AbstractGeometry end
abstract type AbstractDomain end
abstract type AbstractSimplex <: AbstractDomain end
abstract type AbstractCube <: AbstractDomain end
abstract type AbstractPolySet end

import Base: parse
_parse(::Type{<:AbstractString}, x::AbstractString) = x
_parse(::Type{T}, x::AbstractString) where T<:Number = Base.parse(T,x)


include("domain.jl")
include("barycentriccoordinates.jl")
include("referenceelement.jl")
include("symmetryorbits.jl")
include("jacobi.jl")
include("lagrange.jl")
include("monomials.jl")
include("geometry.jl")
include("weights.jl")
include("quadraturerule.jl")
include("test.jl")
include("transform.jl")
include("compactrule.jl")
include("generate.jl")
include("integrate.jl")
include("properties.jl")
include("optimize.jl")
include("dune.jl")

export AbstractCube, AbstractDomain, AbstractGeometry, AbstractSimplex, AffineGeometry,
  AbstractPolySet, BarycentricMonomial, CompactQuadratureRule,
  CompactQuadratureRuleWithWeights, Hexahedron, JacobiPolySet, LagrangePolySet, Line,
  MonomialPolySet, MultiLinearGeometry, Point, Prism, Pyramid, QuadratureRule,
  Quadrilateral, ReferenceElement, SymmetryOrbit, Tetrahedron, Triangle
export ctype, dimension, domain, domaintype, duneReferenceElement, expand, expandall,
  facets, generate, getProperties, getWeights, integrate, isInside, isPositive, isPI,
  optimize, order, position, symmetryOrbits, testQuadratureRule, transform, vertices,
  volume

end # end module FiniteElementQuadratureRules
