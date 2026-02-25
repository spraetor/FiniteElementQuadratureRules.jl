
function testWeights(qr::QuadratureRule)
  check = isapprox(sum(qr.weights), volume(qr.ref), atol=1.e-12)
  if !check
    println("$(sum(qr.weights)) != $(volume(qr.ref))")
  end
  return check
end

function testQuadratureRule(qr::QuadratureRule{Ω,T,P}) where {Ω,T,P}
  polyset = JacobiPolySet(domain(qr),qr.degree)
  max_error = zero(T)
  for (f,I) in zip(polyset.basis, polyset.integrals)
    Q = sum(qr.weights .* f.(qr.points))
    if abs(Q-I) > 1.e-12
      println("$(Q) != $(I)")
    end
    max_error = max(max_error, abs(Q-I))
  end
  return max_error < 1.e-12
end