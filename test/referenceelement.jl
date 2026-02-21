
using FiniteElementQuadratureRules: checkInside, checkStrictlyInside

@testset "ReferenceElement" begin

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
  end

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
