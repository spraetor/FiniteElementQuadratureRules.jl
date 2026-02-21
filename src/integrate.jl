"""
  integrate(f::Function, qr::QuadratureRule)

Compute the integral ∫f(x)dx over a domain `Ω=domain(qr)` using a quadrature rule.
"""
function integrate(f::Function, qr::QuadratureRule)
  sum(f.(qr.points) .* qr.weights)
end