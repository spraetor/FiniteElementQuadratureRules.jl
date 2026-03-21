import YAML

@testset "Compute quadrature weights" begin
  rules_root = joinpath(@__DIR__, "..", "rules")

  @testset "Expand Compact Rule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(data)
    cqrw = CompactQuadratureRuleWithWeights(data)

    qr = expand(cqr)
    qrw = expand(cqrw)

    @test qr.weights ≈ qrw.weights
  end

  @testset "Weights of compact rule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(data)
    qr = expand(cqr)

    w1 = getWeights(cqr.domain, cqr.degree, qr.points, cqr.orbits)
    w2 = getWeights(ReferenceElement(cqr.domain), cqr.degree, qr.points, cqr.orbits)
    @test w1 ≈ w2

    w3 = getWeights(Float64, cqr.domain, cqr.degree, qr.points, cqr.orbits)
    w4 = getWeights(Float64, ReferenceElement(cqr.domain), cqr.degree, qr.points, cqr.orbits)
    @test w1 ≈ w3
    @test w1 ≈ w4

    ow1 = getOrbitWeights(cqr.domain, cqr.degree, qr.points, cqr.orbits)
    ow2 = getOrbitWeights(ReferenceElement(cqr.domain), cqr.degree, qr.points, cqr.orbits)
    ow3 = getOrbitWeights(Float64, cqr.domain, cqr.degree, qr.points, cqr.orbits)
    ow4 = getOrbitWeights(Float64, ReferenceElement(cqr.domain), cqr.degree, qr.points, cqr.orbits)
    @test ow1 ≈ ow2
    @test ow1 ≈ ow3
    @test ow1 ≈ ow4
    @test length(ow1) == sum(cqr.orbits)

    cw = getCompactWeights(cqr)
    @test length(cw) == sum(cqr.orbits)
    @test FiniteElementQuadratureRules.transformWeights(cqr.domain, cw) ≈ ow1
  end

  @testset "Weights of expanded rule" begin
    data = YAML.load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(data)
    qr = expand(cqr)

    w1 = getWeights(domain(qr), qr.degree, qr.points)
    w2 = getWeights(qr.ref, qr.degree, qr.points)
    @test w1 ≈ w2
    @test w1 ≈ qr.weights

    w3 = getWeights(Float64, domain(qr), qr.degree, qr.points)
    w4 = getWeights(Float64, qr.ref, qr.degree, qr.points)
    @test w1 ≈ w3
    @test w1 ≈ w4

    cw = getCompactWeights(qr, cqr.orbits)
    @test FiniteElementQuadratureRules.transformWeights(domain(qr), cw) ≈ getOrbitWeights(qr.ref, qr.degree, qr.points, cqr.orbits)
  end
end
