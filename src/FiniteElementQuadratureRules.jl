module FiniteElementQuadratureRules

abstract type AbstractDomain{D} end

include("domain.jl")
include("generate.jl")
include("quadraturerule.jl")

export AbstractDomain, Point, Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid, dimension
export CompactQuadratureRule, QuadratureRule, generate

end # end module FiniteElementQuadratureRules
