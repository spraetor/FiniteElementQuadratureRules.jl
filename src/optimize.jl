using Optim

# Expand the orbital parameters r in the given symmetry orbit `so`
function _expandorbit(so::SymmetryOrbit, r::AbstractVector)
  @assert length(r) == so.args # number or parameters
  return so.expand(r...)
end

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
      n = so[k].args               # number of orbital parameters
      for _ = 1:orbits[k]          # number of orbits of this type
        r = _clamporbit(so[k],positions[l:l+n-1])
        # Evaluate on the domain reference element, consistent with `getWeights`.
        points = transformcoordinates(domain, collect(_expandorbit(so[k],r)))
        A[i,j] = T(0)
        for x in points
          A[i,j] += pᵢ(x)
        end
        l = l + n
        j = j + 1
      end
    end
  end

  w = A\b
  return sum((A*w - b).^2)
end



function optimize(qr::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain, T<:Real}
  domain = qr.domain
  degree = qr.degree
  orbits = qr.orbits
  positions = qr.positions

  let f = (p) -> _residual(p,(domain,degree,orbits))
    result = Optim.optimize(f, positions; autodiff = :forward)
    println("residuum = $(result.minimum)")
    return CompactQuadratureRule(domain,degree,orbits,result.minimizer)
  end
end
