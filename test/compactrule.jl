using StaticArrays: @SVector

const F = Float64

function test_quadrature_rule(qr::QuadratureRule{Ω,T,P}) where {Ω<:Triangle,T<:Real,P}
  f0 = (x::P) -> 1.0
  f1 = (x::P) -> f0(x) + x[1] + x[2]
  f2 = (x::P) -> f1(x) + x[1]*x[2] + x[1]^2 + x[2]^2
  f3 = (x::P) -> f2(x) + x[1]^2*x[2] + x[1]*x[2]^2 + x[1]^3 + x[2]^3
  f4 = (x::P) -> f3(x) + x[1]^3*x[2] + x[1]^2*x[2]^2 + x[1]*x[2]^3 + x[1]^4 + x[2]^4

  println("$(integrate(f0,qr)), $(integrate(f1,qr)), $(integrate(f2,qr)), $(integrate(f3,qr)), $(integrate(f4,qr))")
end

function pointsAsMatrix(qr::QuadratureRule{Ω,T,P}) where {Ω<:Triangle,T<:Real,P}
  dim = dimension(Ω())
  nPoints = length(qr.points)
  n = length(qr.points[1])

  ref = duneReferenceElement(qr.domain)
  mat = zeros(T,nPoints,dim)
  for i in 1:nPoints
    p = qr.points[i]
    q = @SVector zeros(T,dim)
    for j = 1:n
      q += p[j] * position(ref,j,dim)
    end
    for j in 1:dim
      mat[i,j] = q[j]
    end
  end
  return mat
end

tri = Triangle()

let
  cqr = CompactQuadratureRuleWithWeights(tri, 1, [0,1,0], F[0.0], F[1/3])
  qr = expand(cqr)
  println("QR(1): length=$(length(qr))")
  test_quadrature_rule(qr)
end

let
  cqr = CompactQuadratureRuleWithWeights(tri, 2, [1,2,0], F[0.5,0.0], F[0.45,2/15,0.05])
  qr = expand(cqr)
  println("QR(2): length=$(length(qr))")
  println("  points=$(pointsAsMatrix(qr))")
  test_quadrature_rule(qr)
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
  println("QR(3): length=$(length(qr))")
  println("  points=$(pointsAsMatrix(qr))")
  test_quadrature_rule(qr)
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
  println("QR(4): length=$(length(qr))")
  println("  points=$(pointsAsMatrix(qr))")
  test_quadrature_rule(qr)
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
  println("QR(5): length=$(length(qr))")
  println("  points=$(pointsAsMatrix(qr))")
  test_quadrature_rule(qr)
end