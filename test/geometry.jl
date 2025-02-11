using StaticArrays: SVector

@testset "Geometry" begin

  for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
    ref = ReferenceElement(Domain())

    coordVector = [SVector{length(c)+1,Float64}([c;0]) for c in ref.coordinates]
    geo = MultiLinearGeometry(ref, coordVector)
    println(typeof(geo))

    λ = position(ref,1,0)
    # @test geo(λ) ≈ [λ;0]
  end
end