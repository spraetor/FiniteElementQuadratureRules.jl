
function testWeights(qr::QuadratureRule{Ω,T,P}; tol::Real = sqrt(eps(T))) where {Ω,T,P}
  isapprox(sum(qr.weights), volume(qr.ref), atol=tol)
end

function quadratureAccuracy(qr::QuadratureRule{Ω,T,P}) where {Ω,T,P}
  polyset = JacobiPolySet(domain(qr), qr.degree)
  max_error = zero(T)
  for (f, I) in zip(polyset.basis, polyset.integrals)
    Q = sum(qr.weights .* f.(qr.points))
    max_error = max(max_error, abs(Q - I))
  end
  return max_error
end

function testQuadratureRule(qr::QuadratureRule{Ω,T,P}; tol::Real = sqrt(eps(T))) where {Ω,T,P}
  return quadratureAccuracy(qr) < tol
end
