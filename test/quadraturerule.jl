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

  @testset "Preserve Extra Fields" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    qr = expand(CompactQuadratureRuleWithWeights(data))

    data2 = Dict(qr; extra_fields=Dict("comment" => "keep me"))
    @test data2["comment"] == "keep me"

    mktempdir() do dir
      path = joinpath(dir, "rule.yml")
      write_file(path, qr; reference="TEST", extra_fields=Dict("comment" => "keep me"))
      data3 = YAML.load_file(path)
      @test data3["comment"] == "keep me"
      @test data3["reference"] == "TEST"
    end
  end

  @testset "Quadrature Accuracy" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    qr = expand(CompactQuadratureRuleWithWeights(data))
    tol = sqrt(eps(Float64))
    @test quadratureAccuracy(qr) < tol
    @test testQuadratureRule(qr; tol)
  end

end
