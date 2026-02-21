using FiniteElementQuadratureRules
using Test

@testset "FiniteElementQuadratureRules.jl" begin
    include("domain.jl")
    include("referenceelement.jl")
    include("geometry.jl")
    include("barycentriccoordinates.jl")
    include("transformcoordinates.jl")
    include("polyset.jl")
    include("jacobi.jl")
    include("compactrule.jl")
end
