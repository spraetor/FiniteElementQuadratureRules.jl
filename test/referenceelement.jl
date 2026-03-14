
using FiniteElementQuadratureRules: checkInside, checkStrictlyInside

@testset "ReferenceElement" begin

  @testset "Properties of the reference element" begin
    for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
      domain = Domain()
      ref = ReferenceElement(domain)

      @test domaintype(ref) == Domain
      @test dimension(ref) == dimension(domain)
      @test length(ref.coordinates) == vertices(domain)
      @test length(ref.facets) == facets(domain)

      # Element center (codim 0): inside and, except for Point, strictly inside.
      x0 = position(ref, 1, 0)
      @test checkInside(ref, x0)
      if dimension(domain) > 0
        @test checkStrictlyInside(ref, x0)
      end

      # Test the direct size methods
      @test size(ref, 0) == 1
      @test size(ref, dimension(ref)) == vertices(domain)
      if dimension(ref) > 1
        @test size(ref, 1) == facets(domain)
      end
      if domain == Tetrahedron()
        @test size(ref,2) == 6
      elseif domain == Pyramid()
        @test size(ref,2) == 8
      elseif domain == Prism()
        @test size(ref,2) == 9
      elseif domain == Hexahedron()
        @test size(ref,2) == 12
      end

      # Test the sub-size methods
      @test size(ref, 1,0,0) == 1
      if dimension(ref) > 1
        for i in eachindex(ref.facets)
          @test size(ref, i,1,dimension(ref)) == length(ref.facets[i])
        end
      end

    end
  end

  @testset "CheckInside" begin
    # Centers of boundary sub-entities should lie on the boundary:
    # inside, but not strictly inside.
    for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
      domain = Domain()
      ref = ReferenceElement(domain)
      dim = dimension(domain)

      # check vertices and facets
      codims = (Int[], [1], [1,2], [1,3])[dim+1]

      for c in codims
        n = c == 1 ? facets(domain) : vertices(domain)
        for i in 1:n
          x = position(ref, i, c)
          @test checkInside(ref, x)
          @test !checkStrictlyInside(ref, x)
        end
      end
    end
  end

end
