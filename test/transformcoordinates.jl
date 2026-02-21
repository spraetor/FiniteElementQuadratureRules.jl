using FiniteElementQuadratureRules: transformcoordinates, transformweights

@testset "CoordinateTransforms" begin
  # Line is treated as cube-like 1D reference coordinates.
  @test transformcoordinates(Line(), [[-0.4], [0.2]]) == [[-0.4], [0.2]]

  tri_bary = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.4, 0.2, 0.4],
    [1/3, 1/3, 1/3],
  ]
  tri_ref = transformcoordinates(Triangle(), tri_bary)
  @test tri_ref ≈ [
    [-1.0, -1.0],
    [1.0, -1.0],
    [-1.0, 1.0],
    [-0.6, -0.2],
    [-1/3, -1/3],
  ]

  tet_bary = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
    [0.25, 0.25, 0.25, 0.25],
    [0.35, 0.15, 0.25, 0.25],
  ]
  tet_ref = transformcoordinates(Tetrahedron(), tet_bary)
  @test tet_ref ≈ [
    [-1.0, -1.0, -1.0],
    [1.0, -1.0, -1.0],
    [-1.0, 1.0, -1.0],
    [-1.0, -1.0, 1.0],
    [-0.5, -0.5, -0.5],
    [-0.7, -0.5, -0.5],
  ]

  # Prism: first 3 components are barycentric for triangular base, 4th is z in reference coordinates.
  prism_bary = [[1/3, 1/3, 1/3, 0.25], [1/2, 1/4, 1/4, -0.5]]
  prism_ref = transformcoordinates(Prism(), prism_bary)
  @test prism_ref[1] ≈ [-1/3, -1/3, 0.25]
  @test prism_ref[2] ≈ [-1/2, -1/2, -0.5]

  # Pyramid symmetry orbits are already in 3D reference coordinates.
  py_orbit = first(symmetryOrbits(Float64, Pyramid()))
  py_point = py_orbit.expand(0.2)[1]
  @test transformcoordinates(Pyramid(), [py_point]) == [py_point]

  # Weight transforms.
  @test transformweights(Line(), [1.0, 2.0]) == [1.0, 2.0]
  @test transformweights(Triangle(), [0.25, 0.5]) == [0.5, 1.0]
  @test transformweights(Tetrahedron(), [0.25, 0.5]) == [0.5, 1.0]
  @test transformweights(Quadrilateral(), [0.25, 0.5]) == [0.25, 0.5]
end
