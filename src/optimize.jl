using Optim

"""
    optimize(qr::CompactQuadratureRule)

Optimize the position of the quadrature points given as symmetric orbits in the
compact rule, by minimizing the quadrature residual on a set of polynomial basis
functions. We use basis functions and their integrals given by a `JacobiPolySet`.
"""
function optimize(qr::CompactQuadratureRule)
  domain = qr.domain
  degree = qr.degree
  orbits = qr.orbits
  T = ctype(qr)
  positions = Vector{T}(qr.positions)

  let f = (p) -> _residual(p,(domain,degree,orbits))
    options = Optim.Options(
      x_abstol = sqrt(eps(float(T))),
      g_abstol = sqrt(eps(float(T))),
      f_abstol = eps(float(T)),
      iterations = 10_000,
    )

    # Use a gradient-based method explicitly. The Optim default for scalar
    # objectives is Nelder-Mead, which is too weak here and often stalls early.
    result = Optim.optimize(
      f, positions, BFGS(), options;
      autodiff = :forward,
    )
    return CompactQuadratureRule(domain,degree,orbits,result.minimizer)
  end
end


"""
  optimize(qr::CompactQuadratureRuleWithWeights)

Optimize the position of the quadrature points given as symmetric orbits in the
compact rule, by minimizing the quadrature residual on a set of polynomial basis
functions. We use basis functions and their integrals given by a `JacobiPolySet`.
"""
function optimize(qr::CompactQuadratureRuleWithWeights)
  optimize(CompactQuadratureRule(qr.domain, qr.degree, qr.orbits, qr.positions))
end

# Restrict the arguments given by the vector `r` to the valid range for the
# symmetry orbit `so`.
function _clamporbit(so::SymmetryOrbit, r::AbstractVector)
  return r # TODO
end

# Compute the quadrature residual given by `A_ij * w_j - b_i` with
# - `b_i` the exact integral if the ith basis function pᵢ
# - `w_j` the quadrature weight associated to the jth symmetry orbit
# - `A_ij` the sum ∑ₖ pᵢ(xₖ) for {xₖ} the points in the jth symmetry orbit.
function _residual(positions::AbstractVector{T}, params) where {T}
  domain, degree, orbits = params

  polyset = JacobiPolySet(T, domain, degree)
  so = symmetryOrbits(T,domain)
  nDifferentWeights = sum(orbits)
  A = zeros(T, length(polyset.basis), nDifferentWeights)
  b = zeros(T, length(polyset.basis))

  for i in eachindex(polyset.basis)
    pᵢ = polyset.basis[i]
    b[i] = polyset.integrals[i]

    j = 1
    l = 1
    for k in eachindex(orbits)     # types of symmetry orbits
      n = args(so[k])              # number of orbital parameters
      for _ = 1:orbits[k]          # number of orbits of this type
        r = _clamporbit(so[k],view(positions,l:l+n-1))
        # Evaluate on the domain of the reference element.
        points = transformCoordinates(domain, expand(so[k],r))
        A[i,j] = sum(pᵢ.(points))
        l = l + n
        j = j + 1
      end
    end
  end

  w = A\b
  return sum((A*w - b).^2)
end