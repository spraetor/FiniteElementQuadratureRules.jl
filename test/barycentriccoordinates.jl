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

  let ref = ReferenceElement(Triangle, [[0.0, 0.0], [2.0, 0.0], [0.5, 1.5]], [[1,2], [1,3], [2,3]], 1.5)
    @test barycentricCoordinates(ref, [0.0, 0.0]) ≈ [1.0, 0.0, 0.0]
    @test barycentricCoordinates(ref, [2.0, 0.0]) ≈ [0.0, 1.0, 0.0]
    @test barycentricCoordinates(ref, [0.5, 1.5]) ≈ [0.0, 0.0, 1.0]
    @test barycentricCoordinates(ref, [5/6, 0.5]) ≈ [1/3, 1/3, 1/3]
  end

  let ref = ReferenceElement(Prism,
      [[0.0, 0.0, -2.0], [2.0, 0.0, -2.0], [0.5, 1.5, -2.0],
       [0.0, 0.0, 4.0], [2.0, 0.0, 4.0], [0.5, 1.5, 4.0]],
      [[1,2,4,5], [1,3,4,6], [2,3,5,6], [1,2,3], [4,5,6]], 9.0)
    @test barycentricCoordinates(ref, [0.5, 1.5, -2.0]) ≈ [0.0, 0.0, 1.0, -1.0]
    @test barycentricCoordinates(ref, [5/6, 0.5, 1.0]) ≈ [1/3, 1/3, 1/3, 0.0]
  end
end
