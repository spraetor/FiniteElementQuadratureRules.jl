
@testset "Domain" begin

  @test domaintype(0,"simplex") == Point
  @test domaintype(1,"simplex") == Line
  @test domaintype(2,"simplex") == Triangle
  @test domaintype(2,"cube") == Quadrilateral
  @test domaintype(3,"simplex") == Tetrahedron
  @test domaintype(3,"cube") == Hexahedron
  @test domaintype(3,"prism") == Prism
  @test domaintype(3,"pyramid") == Pyramid

end