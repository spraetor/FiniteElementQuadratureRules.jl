
@testset "Domain" begin

  @test domain(0,"simplex") == Point
  @test domain(1,"simplex") == Line
  @test domain(2,"simplex") == Triangle
  @test domain(2,"cube") == Quadrilateral
  @test domain(3,"simplex") == Tetrahedron
  @test domain(3,"cube") == Hexahedron
  @test domain(3,"prism") == Prism
  @test domain(3,"pyramid") == Pyramid

end