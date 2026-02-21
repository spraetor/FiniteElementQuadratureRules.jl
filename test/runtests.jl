using FiniteElementQuadratureRules
using Test

@testset "FiniteElementQuadratureRules.jl" begin
    include("domain.jl")
    include("referenceelement.jl")
    include("geometry.jl")
    include("polyset.jl")
    include("jacobi.jl")
end

include("compactrule.jl")