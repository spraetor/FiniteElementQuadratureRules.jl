using FiniteElementQuadratureRules: barycentricCoordinates

@testset "BarycentricCoordinates" begin
  # Triangle vertices in reference coordinates.
  let domain = Triangle()
    @test barycentricCoordinates(domain, [-1.0, -1.0]) ≈ [1.0, 0.0, 0.0]
    @test barycentricCoordinates(domain, [1.0, -1.0]) ≈ [0.0, 1.0, 0.0]
    @test barycentricCoordinates(domain, [-1.0, 1.0]) ≈ [0.0, 0.0, 1.0]
  end

  # Tetrahedron vertices in reference coordinates.
  let domain = Tetrahedron()
    @test barycentricCoordinates(domain, [-1.0, -1.0, -1.0]) ≈ [1.0, 0.0, 0.0, 0.0]
    @test barycentricCoordinates(domain, [1.0, -1.0, -1.0]) ≈ [0.0, 1.0, 0.0, 0.0]
    @test barycentricCoordinates(domain, [-1.0, 1.0, -1.0]) ≈ [0.0, 0.0, 1.0, 0.0]
    @test barycentricCoordinates(domain, [-1.0, -1.0, 1.0]) ≈ [0.0, 0.0, 0.0, 1.0]
  end

  # Prism: triangular base in barycentric coordinates + reference z coordinate.
  let domain = Prism()
    @test barycentricCoordinates(domain, [-1.0, -1.0, -0.5]) ≈ [1.0, 0.0, 0.0, -0.5]
    @test barycentricCoordinates(domain, [1.0, -1.0, 0.25]) ≈ [0.0, 1.0, 0.0, 0.25]
    @test barycentricCoordinates(domain, [-1.0, 1.0, 0.75]) ≈ [0.0, 0.0, 1.0, 0.75]
    @test barycentricCoordinates(domain, [-1/3, -1/3, 0.1]) ≈ [1/3, 1/3, 1/3, 0.1]
  end

  for domain in (Point(),Line(),Quadrilateral(),Hexahedron(),Pyramid())
    ref = ReferenceElement(domain)
    x = position(ref,1,0)
    @test barycentricCoordinates(domain, x) == x
  end
end
