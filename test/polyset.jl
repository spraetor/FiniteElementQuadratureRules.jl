using FiniteElementQuadratureRules: _allexponents

@testset "MonomialPolySet" begin

  # check exponent tuples
  for len in 1:4
    @test length(_allexponents(len,0)) == 1
    @test length(_allexponents(len,1)) == len

    for total in 0:5
      tuples = _allexponents(len,total)
      @test length(tuples) == binomial(total+len-1, len-1)
      @test all(length.(tuples) .== len)
      @test all(sum.(tuples) .== total)
    end
  end

  tri = Triangle()
  p = MonomialPolySet(tri, 5)
  @test integrate(p.basis[1], tri) ≈ 1.0
end