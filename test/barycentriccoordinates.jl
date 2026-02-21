using FiniteElementQuadratureRules: barycentriccoordinates

@testset "BarycentricCoordinates" begin
  # Triangle vertices in reference coordinates.
  let domain = Triangle()
    @test barycentriccoordinates(domain, [-1.0, -1.0]) ≈ [1.0, 0.0, 0.0]
    @test barycentriccoordinates(domain, [1.0, -1.0]) ≈ [0.0, 1.0, 0.0]
    @test barycentriccoordinates(domain, [-1.0, 1.0]) ≈ [0.0, 0.0, 1.0]
  end

  # Tetrahedron vertices in reference coordinates.
  let domain = Tetrahedron()
    @test barycentriccoordinates(domain, [-1.0, -1.0, -1.0]) ≈ [1.0, 0.0, 0.0, 0.0]
    @test barycentriccoordinates(domain, [1.0, -1.0, -1.0]) ≈ [0.0, 1.0, 0.0, 0.0]
    @test barycentriccoordinates(domain, [-1.0, 1.0, -1.0]) ≈ [0.0, 0.0, 1.0, 0.0]
    @test barycentriccoordinates(domain, [-1.0, -1.0, 1.0]) ≈ [0.0, 0.0, 0.0, 1.0]
  end

  # Prism: triangular base in barycentric coordinates + reference z coordinate.
  let domain = Prism()
    @test barycentriccoordinates(domain, [-1.0, -1.0, -0.5]) ≈ [1.0, 0.0, 0.0, -0.5]
    @test barycentriccoordinates(domain, [1.0, -1.0, 0.25]) ≈ [0.0, 1.0, 0.0, 0.25]
    @test barycentriccoordinates(domain, [-1.0, 1.0, 0.75]) ≈ [0.0, 0.0, 1.0, 0.75]
    @test barycentriccoordinates(domain, [-1/3, -1/3, 0.1]) ≈ [1/3, 1/3, 1/3, 0.1]
  end
end
