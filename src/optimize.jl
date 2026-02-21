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
        points = transformCoordinates(domain, SVector(expand(so[k],r)))
        A[i,j] = sum(pᵢ.(points))
        l = l + n
        j = j + 1
      end
    end
  end

  w = A\b
  return sqrt(sum((A*w - b).^2))
end



function optimize(qr::CompactQuadratureRule)
  domain = qr.domain
  degree = qr.degree
  orbits = qr.orbits
  positions = qr.positions

  let f = (p) -> _residual(p,(domain,degree,orbits))
    result = Optim.optimize(f, positions; autodiff = :forward)
    return CompactQuadratureRule(domain,degree,orbits,result.minimizer)
  end
end
