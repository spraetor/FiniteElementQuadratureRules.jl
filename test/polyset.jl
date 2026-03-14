using FiniteElementQuadratureRules: _allexponents

function test_polyset(ps::AbstractPolySet)
  @test length(ps.basis) == length(ps.integrals)
  ref = ReferenceElement(ps.domain)
  for f in ps.basis
    f(position(ref,1,0)) # evaluate in the center of the domain
  end
end

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
  # test_polyset(p)
  @test integrate(p.basis[1], tri) ≈ 1.0
end

@testset "LagrangePolySet" begin

  for degree in 0:3
    lps2 = LagrangePolySet(Triangle(), degree)
    test_polyset(lps2)

    lps3 = LagrangePolySet(Tetrahedron(), degree)
    test_polyset(lps3)
  end

end