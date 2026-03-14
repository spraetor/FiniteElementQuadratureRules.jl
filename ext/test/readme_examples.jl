import YAML

@testset "README Examples" begin
  rules_root = joinpath(@__DIR__, "..", "..", "rules")

  @testset "Transform Rule To Dune Convention" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    qr = expand(CompactQuadratureRule(Float64, data))

    ref_dune = duneReferenceElement(domain(qr))
    qr_dune = transform(qr, ref_dune)
    @test qr_dune.degree == qr.degree
    @test length(qr_dune) == length(qr)
  end
end
