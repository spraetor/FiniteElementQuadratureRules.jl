using FiniteElementQuadratureRules: _jacobi

@testset "Jacobi" begin
  for n in 0:5
    for a = 0:5
      for b in 0:5
        @test _jacobi(n,a,b,1.0) == binomial(n+a,n)

        X = range(-1,1,7)
        @test _jacobi(n,a,b,X) ≈ map((x) -> _jacobi(n,a,b,x), X)
        @test _jacobi(n,a,b,-X) ≈ (-1)^n .* _jacobi(n,b,a,X)

        if n > 0 && b > 0
          @test (2n+a+b).*_jacobi(n,a,b-1,X) ≈ (n+a+b).*_jacobi(n,a,b,X) .+ (n+a).*_jacobi(n-1,a,b,X)
        end
      end
    end
  end

  # 5-point Gauss-Legendre rule on [-1,1].
  gl_x = (-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640)
  gl_w = (0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891)
  integrate_line(f) = sum(w * f((x,)) for (x,w) in zip(gl_x, gl_w))

  for degree in 0:8
    jac = JacobiPolySet(Line(), degree)
    @test jac.integrals ≈ map(f -> integrate_line(f), jac.basis) atol=1e-13
  end

  let tri = Triangle()
    cqr = CompactQuadratureRuleWithWeights(tri, 5, [0,4,3],
      [ 0.00000000000000000000000000000000000,
        0.05752768441141010566081751776543190,
        0.25685910726195907606389083148220168,
        0.45783683807916110193850321742406830,
        0.00000000000000000000000000000000000,
        0.36329807415368604570550633618418105,
        0.00000000000000000000000000000000000,
        0.13226458163271398535388822004364736,
        0.22100121875989000797812820146484192,
        0.07819258362551702199888597846982583 ],
      [ 0.00141884794135849195859201314216757,
        0.02325227091923514227899684743112030,
        0.09180247526152571475403829582713371,
        0.06906086075456558705677702806082523,
        0.01238113000735325822823625385145974,
        0.00696115728097842131688535505330661,
        0.05455715193999251909734296553127691 ])
    qr = expand(cqr)

    for degree in 0:5
      jac = JacobiPolySet(tri, degree)
      @test jac.integrals ≈ map(f -> integrate(f, qr), jac.basis) atol=1e-12
    end
  end

  let tet = Tetrahedron()
    cqr = CompactQuadratureRule(tet, 14, [0, 5, 0, 4, 7],
    [ 0.3272533625238485639093096692685289,
      0.0447613044666850808837942096478842,
      0.0861403311024363536537208740298857,
      0.2087626425004322968265357083976176,
      0.0141049738029209600635879152102928,
      0.1021653241807768123476692526982584,
      0.5739463675943338202814002893460107,
      0.4075700516600107157213295651301783,
      0.0922278701390201300000000000000000,
      0.0156640007402803585557586709578084,
      0.7012810959589440327139967673208426,
      0.2254963562525029053780724154201103,
      0.4769063974420887115860583354107011,
      0.3905984281281458000000000000000000,
      0.2013590544123922168123077327235092,
      0.0161122880710300298578026931548371,
      0.1061350679989021455556139029848079,
      0.0327358186817269284944004077912660,
      0.0035979076537271666907971523385925,
      0.5636383731697743896896816630648502,
      0.2302920722300657454502526874135652,
      0.1907199341743551862712487790637898,
      0.3676255095325860844092206775991167,
      0.2078851380230044950717102125250735,
      0.3312104885193449000000000000000000,
      0.7192323689817295295023401840796991,
      0.1763279118019329762157993033636973,
      0.0207602362571310090754973440611644,
      0.5278249952152987298409240075817276,
      0.4372890892203418165526238760841918,
      0.0092201651856641949463177554949220,
      0.5483674544948190728994910505607746,
      0.3447815506171641228703671870920331,
      0.0867217283322215394629438740085828 ])
    qr = expand(cqr)
    for degree in 0:4
      jac = JacobiPolySet(tet, degree)
      @test jac.integrals ≈ map(f -> integrate(f, qr), jac.basis) atol=1e-11
    end
  end

  # Cube domains: tensor-product Gauss integration.
  function integrate_quad(f)
    s = 0.0
    for (x,wx) in zip(gl_x, gl_w), (y,wy) in zip(gl_x, gl_w)
      s += wx*wy*f((x,y))
    end
    s
  end
  function integrate_hex(f)
    s = 0.0
    for (x,wx) in zip(gl_x, gl_w), (y,wy) in zip(gl_x, gl_w), (z,wz) in zip(gl_x, gl_w)
      s += wx*wy*wz*f((x,y,z))
    end
    s
  end

  for degree in 0:6
    jac = JacobiPolySet(Quadrilateral(), degree)
    @test jac.integrals ≈ map(f -> integrate_quad(f), jac.basis) atol=1e-12
  end
  for degree in 0:6
    jac = JacobiPolySet(Hexahedron(), degree)
    @test jac.integrals ≈ map(f -> integrate_hex(f), jac.basis) atol=1e-11
  end

  # Prism and pyramid currently lack robust end-to-end rule support in tests;
  # at least verify the analytic constants stored in the Jacobi polyset.
  @test JacobiPolySet(Prism(), 0).integrals[1] ≈ 4*sqrt(2.0)
  @test JacobiPolySet(Pyramid(), 0).integrals[1] ≈ 16/3
end
