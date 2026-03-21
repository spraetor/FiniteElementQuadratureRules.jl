import YAML

@testset "QuadratureRule" begin
  rules_root = joinpath(@__DIR__, "..", "rules")

  @testset "Construct QuadratureRuleRule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRuleWithWeights(data)
    qr = expand(cqr)

    qr1 = QuadratureRule(domain(qr), qr.degree, qr.points)
    qr2 = QuadratureRule(domain(qr), qr.degree, qr.points, qr.weights)

    @test qr.weights ≈ qr1.weights
    @test qr.weights ≈ qr2.weights
    @test qr.points ≈ qr1.points
    @test qr.points ≈ qr2.points
  end

  @testset "Export QuadratureRuleRule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRuleWithWeights(data)
    qr = expand(cqr)

    data2 = Dict(qr)
    @test data["degree"] == data2["degree"]
    @test data["dim"] == data2["dim"]
    @test data["region"] == data2["region"]
    @test haskey(data2, "accuracy")
    @test parse(Float64, data2["accuracy"]) ≈ quadratureAccuracy(qr)
  end

  @testset "Quadrature Accuracy" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    qr = expand(CompactQuadratureRuleWithWeights(data))
    tol = sqrt(eps(Float64))
    @test quadratureAccuracy(qr) < tol
    @test testQuadratureRule(qr; tol)
  end

  @testset "Read expanded QuadratureRuleRule" begin
    data = YAML.load_file(joinpath(rules_root, "expanded", "WV15", "tetrahedron", "1-1.yml"))
    qr = QuadratureRule(data)
    # @test testQuadratureRule(qr)

    qr2 = QuadratureRule(Float64,data)
    # @test testQuadratureRule(qr2)
  end
end
