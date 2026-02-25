import NonlinearSolve as NLS

"""
    optimize(qr::CompactQuadratureRule)

Optimize the position of the quadrature points given as symmetric orbits in the
compact rule, by minimizing the quadrature residual on a set of polynomial basis
functions. We use basis functions and their integrals given by a `JacobiPolySet`.
"""
function optimize(qr::CompactQuadratureRule{Ω,T};
                  maxiters=1_000, show_trace=false) where {Ω,T}
  domain = qr.domain
  orbits = qr.orbits
  polyset = JacobiPolySet(T, domain, qr.degree)
  b = polyset.integrals

  let prob = NLS.NonlinearProblem(_residual, qr.positions, (domain,orbits,polyset,b))
    if show_trace
      # Work around a LinearSolve QR cache type mismatch for non-BLAS element types
      # (e.g. Float128, Double64) by using Julia's native `\` linsolve path.
      linsolve = T <: Union{Float32, Float64} ? nothing : (\)
      sol = NLS.solve(prob, NLS.LevenbergMarquardt(; linsolve),
        show_trace=Val(true),
        maxiters=maxiters)
    else
      sol = NLS.solve(prob, NLS.SimpleTrustRegion(),
        maxiters=maxiters)
    end
    return CompactQuadratureRule(domain,qr.degree,orbits,sol.u)
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
  domain, orbits, polyset, b = params
  so = symmetryOrbits(T,domain)
  nDifferentWeights = sum(orbits)
  A = zeros(T, length(polyset.basis), nDifferentWeights)

  for i in eachindex(polyset.basis)
    pᵢ = polyset.basis[i]

    j = 1
    l = 1
    for k in eachindex(orbits)     # types of symmetry orbits
      n = args(so[k])              # number of orbital parameters
      for _ = 1:orbits[k]          # number of orbits of this type
        r = view(positions,l:l+n-1)
        # Evaluate on the domain of the reference element.
        points = transformCoordinates(domain, expand(so[k],r))
        A[i,j] = sum(pᵢ.(points))
        l = l + n
        j = j + 1
      end
    end
  end

  w = A\b
  return A*w - b
end
