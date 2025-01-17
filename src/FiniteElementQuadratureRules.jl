module FiniteElementQuadratureRules

abstract type AbstractDomain{D} end

include("domain.jl")
include("quadraturerule.jl")

include("generate.jl")

export AbstractDomain, Point, Line, Triangle, Quadrilateral, Tetrahedron, Hexahedron, Prism, Pyramid, dimension
export CompactQuadratureRule, QuadratureRule, generate

end # end module FiniteElementQuadratureRules
