module FiniteElementQuadratureRules

abstract type AbstractLocalBasis end
abstract type AbstractGeometry end
abstract type AbstractDomain end
abstract type AbstractSimplex <: AbstractDomain end
abstract type AbstractCube <: AbstractDomain end

include("domain.jl")
include("referenceelement.jl")
include("lagrange.jl")
include("geometry.jl")
include("quadraturerule.jl")
include("generate.jl")
include("transform.jl")

export AbstractCube, AbstractDomain, AbstractGeometry, AbstractLocalBasis,
  AbstractSimplex, AffineGeometry, CompactQuadratureRule, Hexahedron, Lagrange,
  Line, MultiLinearGeometry, Point, Prism, Pyramid, QuadratureRule, Quadrilateral,
  Tetrahedron, Triangle
export ctype, dimension, domain, domaintype, expand, facets, generate, order,
  transform, vertices, volume

end # end module FiniteElementQuadratureRules
