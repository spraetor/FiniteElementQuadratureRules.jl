module FiniteElementQuadratureRules

abstract type AbstractDomain end
abstract type AbstractSimplex <: AbstractDomain end
abstract type AbstractCube <: AbstractDomain end

include("domain.jl")
include("quadraturerule.jl")
include("generate.jl")

export AbstractDomain, Point, Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid, dimension
export CompactQuadratureRule, QuadratureRule, generate

end # end module FiniteElementQuadratureRules
