module FiniteElementQuadratureRules

abstract type AbstractGeometry end
abstract type AbstractDomain end
abstract type AbstractSimplex <: AbstractDomain end
abstract type AbstractCube <: AbstractDomain end
abstract type AbstractPolySet end

import Base: parse
_parse(::Type{<:AbstractString}, x::AbstractString) = x
_parse(::Type{T}, x::AbstractString) where T<:Number = Base.parse(T,x)

# Function which are available in the FiniteElementQuadratureRulesExportExt
# extension, when BibFormatter, BibInternal, BibParser, and OteraEngine are
# loaded.
function default_chooser end
function duneReferenceElement end
function generate end
function generate_rule_overview end


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

export AbstractCube, AbstractDomain, AbstractGeometry, AbstractSimplex, AffineGeometry,
  AbstractPolySet, BarycentricMonomials, CompactQuadratureRule,
  CompactQuadratureRuleWithWeights, Hexahedron, JacobiPolySet, LagrangePolySet, Line,
  MonomialPolySet, MultiLinearGeometry, Point, Prism, Pyramid, QuadratureRule,
  Quadrilateral, ReferenceElement, SymmetryOrbit, Tetrahedron, Triangle
export args, compact, ctype, default_chooser, dimension, domain, domaintype, duneReferenceElement, expand, expandall,
  facets, generate, generate_rule_overview, getProperties, getQuality, getWeights, integrate, isInside, isPositive, isPI,
  optimize, position, region, symmetryOrbits, testQuadratureRule, testWeights, transform,
  vertices, volume, write_file

end # end module FiniteElementQuadratureRules
