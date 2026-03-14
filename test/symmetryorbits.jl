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
