using Optim
using StaticArrays: SVector

function _clamporbit(so::SymmetryOrbit, r::AbstractVector)
  return r # TODO
end

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
        # Evaluate on the domain reference element, consistent with `getWeights`.
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
