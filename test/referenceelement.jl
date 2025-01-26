
@testset "ReferenceElement" begin

  for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
    domain = Domain()
    ref = ReferenceElement(domain)

    @test domaintype(ref) == Domain
    @test dimension(ref) == dimension(domain)
    @test length(ref.coordinates) == vertices(domain)
    @test length(ref.facets) == facets(domain)
  end

end
