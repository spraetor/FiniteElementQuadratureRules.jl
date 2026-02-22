using StaticArrays: @SVector

const F = Float64

function test_quadrature_rule(qr::QuadratureRule{Ω,T,P}) where {Ω,T,P}
  polyset = JacobiPolySet(qr.domain,qr.degree)
  @test sum(qr.weights) ≈ volume(ReferenceElement(qr.domain))
  max_error = zero(T)
  for (f,I) in zip(polyset.basis, polyset.integrals)
    Q = sum(qr.weights .* f.(qr.points))
    @test Q ≈ I atol=sqrt(eps(T))
    max_error = max(max_error, abs(Q-I))
  end
  return max_error
end

tri = Triangle()
line = Line()

let
  cqr = CompactQuadratureRuleWithWeights(line, 1, [1,0], F[], F[2.0])
  qr = expand(cqr)

  cqr2 = CompactQuadratureRule(line, 1, [1,0], F[])
  qr2 = expand(cqr2)

  @test qr.weights ≈ qr2.weights
  test_quadrature_rule(qr)
  test_quadrature_rule(qr2)
end

let
  cqr = CompactQuadratureRuleWithWeights(tri, 1, [0,1,0], F[0.0], F[1/3])
  qr = expand(cqr)

  cqr2 = CompactQuadratureRule(tri, 1, [0,1,0], F[0.0])
  qr2 = expand(cqr2)

  @test qr.weights ≈ qr2.weights
  test_quadrature_rule(qr)
end

let
  cqr = CompactQuadratureRuleWithWeights(tri, 2, [1,2,0], F[0.5,0.0], F[0.45,2/15,0.05])
  qr = expand(cqr)

  cqr2 = CompactQuadratureRule(tri, 2, [1,2,0], F[0.5,0.0])
  qr2 = expand(cqr2)

  test_quadrature_rule(qr)
  test_quadrature_rule(qr2)
  # @test qr.weights ≈ qr2.weights
end

let
  cqr = CompactQuadratureRuleWithWeights(tri, 3, [0,2,1],
    F[ 0.00000000000000000000000000000000000,
      0.20734517566359092426182782125527331,
      0.00000000000000000000000000000000000,
      0.70653044409095980961019599556083747 ],
    F[ 0.01487291302482058169305106735100205,
      0.22077705784041073851802413812368105,
      0.04884168123405100656112906392932512 ])
  qr = expand(cqr)

  cqr2 = CompactQuadratureRule(tri, 3, [0,2,1],
    F[ 0.00000000000000000000000000000000000,
      0.20734517566359092426182782125527331,
      0.00000000000000000000000000000000000,
      0.70653044409095980961019599556083747 ])
  qr2 = expand(cqr2)

  test_quadrature_rule(qr)
  test_quadrature_rule(qr2)
  # @test qr.weights ≈ qr2.weights
end

let
  cqr = CompactQuadratureRuleWithWeights(tri, 4, [0,4,1],
    F[ 0.50000000000000000000000000000000000,
      0.00000000000000000000000000000000000,
      0.42476396172581058836120087520218113,
      0.13079159382974496719435468035337442,
      0.00000000000000000000000000000000000,
      0.21132486540518711774542560974902127 ],
    F[ 0.02539682539682539682539682539682540,
      0.00634920634920634920634920634920635,
      0.15756242893878361853822583683016899,
      0.10116772979137511162050432189998974,
      0.02142857142857142857142857142857143 ])
  qr = expand(cqr)

  cqr2 = CompactQuadratureRule(tri, 4, [0,4,1],
  F[ 0.50000000000000000000000000000000000,
    0.00000000000000000000000000000000000,
    0.42476396172581058836120087520218113,
    0.13079159382974496719435468035337442,
    0.00000000000000000000000000000000000,
    0.21132486540518711774542560974902127 ])
  qr2 = expand(cqr2)

  # @test qr.weights ≈ qr2.weights
  test_quadrature_rule(qr)
  test_quadrature_rule(qr2)
end

let
  cqr = CompactQuadratureRuleWithWeights(tri, 5, [0,4,3],
    F[ 0.00000000000000000000000000000000000,
      0.05752768441141010566081751776543190,
      0.25685910726195907606389083148220168,
      0.45783683807916110193850321742406830,
      0.00000000000000000000000000000000000,
      0.36329807415368604570550633618418105,
      0.00000000000000000000000000000000000,
      0.13226458163271398535388822004364736,
      0.22100121875989000797812820146484192,
      0.07819258362551702199888597846982583 ],
    F[ 0.00141884794135849195859201314216757,
      0.02325227091923514227899684743112030,
      0.09180247526152571475403829582713371,
      0.06906086075456558705677702806082523,
      0.01238113000735325822823625385145974,
      0.00696115728097842131688535505330661,
      0.05455715193999251909734296553127691 ])
  qr = expand(cqr)


  cqr2 = CompactQuadratureRule(tri, 5, [0,4,3],
  F[ 0.00000000000000000000000000000000000,
    0.05752768441141010566081751776543190,
    0.25685910726195907606389083148220168,
    0.45783683807916110193850321742406830,
    0.00000000000000000000000000000000000,
    0.36329807415368604570550633618418105,
    0.00000000000000000000000000000000000,
    0.13226458163271398535388822004364736,
    0.22100121875989000797812820146484192,
    0.07819258362551702199888597846982583 ])
  qr2 = expand(cqr2)

  # @test qr.weights ≈ qr2.weights
  test_quadrature_rule(qr)
  test_quadrature_rule(qr2)
end

let
  cqr = CompactQuadratureRule(Tetrahedron(), 14, [0, 5, 0, 4, 7],
  F[ 0.3272533625238485639093096692685289,
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
  @test length(qr) == 236
  test_quadrature_rule(qr)
end

let
  # Exactly representable orbit parameters on triangle should be optimizer fixed points.
  cqr = CompactQuadratureRule(tri, 1, [0,1,0], F[0.0])
  oqr = optimize(cqr)
  @test cqr.positions ≈ oqr.positions atol=1e-12
end

let
  # Another minimal case with binary-representable values.
  cqr = CompactQuadratureRule(tri, 2, [1,2,0], F[0.5,0.0])
  oqr = optimize(cqr)
  @test cqr.positions ≈ oqr.positions atol=1e-12
end

let F = Float64
  # Test the optimization on a more advanced example
  cqr = CompactQuadratureRule(tri, 5, [1,2,0],
  F[ 4.7014206410e-01,
     1.0128650732e-01])
  qr = expand(cqr)

  ocqr = optimize(cqr)
  oqr = expand(ocqr)

  @test getProperties(qr) == getProperties(oqr)
  @test cqr.positions ≈ ocqr.positions atol=1e-9

  test_quadrature_rule(oqr)
end


let domain = Quadrilateral()
  # Test compact rules on a Quadrilateral
  cqr = CompactQuadratureRule(domain, 5, [0,1,1],
    F[ 6.83130051063973225548e-01,
       8.81917103688196863500e-01 ])
  qr = expand(cqr)
  oqr = optimize(cqr)

  test_quadrature_rule(qr)
end

let domain = Tetrahedron()
  # Test compact rules on a Tetrahedron
  cqr = CompactQuadratureRule(domain, 5, [0,2,1,0],
    F[ 3.10885919263300609797e-01,
       9.27352503108912264023e-02,
       4.55037041256496494918e-02 ])
  qr = expand(cqr)
  oqr = optimize(cqr)

  test_quadrature_rule(qr)
end

let domain = Pyramid()
  # Test compact rules on a Pyramid
  cqr = CompactQuadratureRule(domain, 5, [3,1,2,0],
    F[ 4.59715761565013385861e-01,
      -3.99197958372461985933e-01,
      -9.99999987016455692414e-01,
       7.06526031546324574207e-01,
      -7.50000000000000000000e-01,
       7.05117122778827601810e-01,
      -8.77776185875954071084e-01,
       4.32882864103540976850e-01,
      -1.52797325760550388420e-01 ])
  qr = expand(cqr)
  oqr = optimize(cqr)

  test_quadrature_rule(qr)
end

let domain = Prism()
  # Test compact rules on a Prism
  cqr = CompactQuadratureRule(domain, 5, [1,0,1,2,0,0],
    F[ 4.87299964550245666116e-01,
       4.45598104670372003762e-01,
       8.71002934865444052729e-01,
       1.00858945982708530915e-01,
       5.70426980705159272061e-01 ])
  qr = expand(cqr)
  oqr = optimize(cqr)

  test_quadrature_rule(qr)
end

let domain = Hexahedron()
  # Test compact rules on a Hexahedron
  cqr = CompactQuadratureRule(domain, 5, [0,1,1,0,0,0,0],
    F[ 7.95822425754221463264e-01,
       7.58786910639328146269e-01 ])
  qr = expand(cqr)
  oqr = optimize(cqr)

  test_quadrature_rule(qr)
end