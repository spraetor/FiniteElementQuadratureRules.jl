import YAML

@testset "README Examples" begin
  rules_root = joinpath(@__DIR__, "..", "rules")

  @testset "Expand Compact Rule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(BigFloat, data)
    qr = expand(cqr)

    @test qr.degree == 4
    @test length(qr) == 6

    mktempdir() do tmp
      out = joinpath(tmp, "4-6.yml")
      write_file(out, cqr; reference=data["reference"])
      write_file(out, qr; reference=data["reference"])
      YAML.write_file(out, Dict(qr; reference=data["reference"]))
      @test isfile(out)
      exported = YAML.load_file(out)
      @test Int(exported["degree"]) == 4
      @test length(exported["coordinates"]) == 6
    end
  end

  @testset "Optimize Rule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(BigFloat, data)
    oqr = optimize(cqr)
    @test oqr.degree == cqr.degree
    @test length(oqr.positions) == length(cqr.positions)
  end
end
