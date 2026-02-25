using YAML: load_file

@testset "Dune reference element" begin
  rules_root = joinpath(@__DIR__, "..", "rules")

  data = load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
  cqr = CompactQuadratureRule(Float64, data)
  oqr = optimize(cqr)
  qr = expand(oqr)
  @test testQuadratureRule(qr)

  dune_qr = transform(qr, duneReferenceElement(domain(qr)))
  @test testWeights(dune_qr)
  # @test testQuadratureRule(dune_qr)
end