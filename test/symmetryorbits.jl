using StaticArrays: SVector

@testset "Triangle symmetry orbits" begin
  orbits = symmetryOrbits(Float64, Triangle())

  @test length(orbits) == 8
  @test args.(orbits) == [0, 1, 2, 2, 2, 2, 2, 3]
  @test length.(orbits) == [1, 3, 6, 3, 2, 2, 2, 1]

  p = SVector(0.2, 0.3)
  c = 1.0 - sum(p)

  @test expand(orbits[5], p) == SVector(
    SVector(0.2, 0.3, c),
    SVector(0.3, 0.2, c),
  )
  @test expand(orbits[6], p) == SVector(
    SVector(0.2, 0.3, c),
    SVector(c, 0.3, 0.2),
  )
  @test expand(orbits[7], p) == SVector(
    SVector(0.2, 0.3, c),
    SVector(0.2, c, 0.3),
  )
end

@testset "Orbit detection" begin
  domain = Triangle()
  points = [
    SVector(0.1, 0.2, 0.7),
    SVector(0.1, 0.7, 0.2),
    SVector(0.4, 0.3, 0.3),
  ]
  weights = [0.125, 0.125, 0.2]

  detected = detectSymmetryOrbits(domain, points;
    weights=weights,
    candidate_orbits=[7, 8])

  @test detected.orbits == [0, 0, 0, 0, 0, 0, 1, 1]
  @test detected.positions ≈ [0.1, 0.2, 0.4, 0.3, 0.3]
  @test detected.weights ≈ [0.125, 0.2]

  asymmetric = detectSymmetryOrbits(domain, points;
    candidate_orbits=[8])
  @test asymmetric.orbits == [0, 0, 0, 0, 0, 0, 0, 3]
  @test length(asymmetric.positions) == 9
end

@testset "Hexahedron identity orbit" begin
  orbits = symmetryOrbits(Float64, Hexahedron())

  @test length(orbits) == 8
  @test args.(orbits) == [0, 1, 1, 1, 2, 2, 3, 3]
  @test length.(orbits) == [1, 6, 8, 12, 24, 24, 48, 1]

  p = SVector(0.2, -0.3, 0.4)
  @test expand(orbits[8], p) == [SVector(0.2, -0.3, 0.4)]
  @test compact(orbits[8], p) == (0.2, -0.3, 0.4)
end
