using YAML
using Base: Filesystem

@testset "Compact rules" begin
  rules_root = joinpath(@__DIR__, "..", "rules", "compact")

  for (root, _, files) in Filesystem.walkdir(rules_root)
    for file in (f for f in files if endswith(f, ".yml"))
      occursin(".opt.yml", file) && continue
      data = YAML.load_file(joinpath(root, file))
      qr = if haskey(data, "coordinates")
        QuadratureRule(Float64, data)
      elseif haskey(data, "weights")
        expand(CompactQuadratureRuleWithWeights(Float64, data))
      else
        expand(CompactQuadratureRule(Float64, data))
      end

      tol = sqrt(eps(Float64))
      @test testWeights(qr; tol)
      @test testQuadratureRule(qr; tol)
    end
  end
end
