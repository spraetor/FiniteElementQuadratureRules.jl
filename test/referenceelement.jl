
@testset "ReferenceElement" begin

  for Domains in Base.uniontypes(AllDomains)
    @test typeof(ReferenceElement(Domain).domain) == Domain
  end

end
