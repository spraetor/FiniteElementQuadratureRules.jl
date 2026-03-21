using YAML: load_file

include(joinpath(@__DIR__, "..", "rules", "data", "taylor_wingate_bos_2007", "convert.jl"))

@testset "TWB05 conversion" begin
  compact_dir = joinpath(@__DIR__, "..", "rules", "compact", "TWB05", "triangle")
  @test isdir(compact_dir)

  source_rules = load_twb05_source_rules(joinpath(@__DIR__, "..", "rules", "data", "taylor_wingate_bos_2007"))
  compact_files = sort(filter(f -> occursin(r"/\d+-\d+\.yml$", f), readdir(compact_dir; join=true)))

  @test length(source_rules) == length(compact_files)

  source_by_degree = Dict(rule.degree => rule for rule in source_rules)
  optimized_degrees = Set([7, 14])
  for file in compact_files
    data = load_file(file)
    @test data["reference"] == "TWB05"
    cqr = CompactQuadratureRuleWithWeights(BigFloat, data)
    qr = expand(cqr)
    @test testWeights(qr; tol = sqrt(eps(Float64)))
    @test testQuadratureRule(qr; tol = sqrt(eps(Float64)))
    if cqr.degree in optimized_degrees
      @test !compare_twb05_rule(source_by_degree[cqr.degree], cqr; tol = 1e-8)
    else
      @test compare_twb05_rule(source_by_degree[cqr.degree], cqr; tol = 1e-8)
    end
  end
end
