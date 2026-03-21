using FiniteElementQuadratureRules
using Test

@testset "FiniteElementQuadratureRules.jl" begin
    include("barycentriccoordinates.jl")
    include("compactrule.jl")
    include("domain.jl")
    include("geometry.jl")
    include("jacobi.jl")
    include("polyset.jl")
    include("quadraturerule.jl")
    include("readme_examples.jl")
    include("referenceelement.jl")
    include("rules.jl")
    include("symmetryorbits.jl")
    include("transformcoordinates.jl")
    include("weights.jl")
end

# run the test in the extension
include(joinpath(@__DIR__, "..", "ext", "test", "runtest.jl"))
