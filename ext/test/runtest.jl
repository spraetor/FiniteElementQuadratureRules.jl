using FiniteElementQuadratureRules
using Test
include(joinpath(@__DIR__, "..", "dune.jl"))
include(joinpath(@__DIR__, "..", "generate_common.jl"))

@testset "FiniteElementQuadratureRulesExportExt.jl" begin
    include("common.jl")
    include("dune.jl")
    include("readme_examples.jl")
end
