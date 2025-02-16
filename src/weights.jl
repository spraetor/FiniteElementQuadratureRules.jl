function getWeights(::Type{T}, domain::Ω, degree::Integer, points::AbstractVector{P}) where {T<:Real,Ω<:AbstractDomain,P<:AbstractVector}

  polyset = PolySet(domain, degree)
  A = zeros(T, length(polyset.basis), length(points))
  b = zeros(T, length(polyset.basis))

  ref = ReferenceElement(domain)
  for i in eachindex(polyset.basis)
    f = polyset.basis[i]
    for j in eachindex(points)
      A[i,j] = f(points[j])
    end
    b[i] = T(volume(ref) * integrate(f,domain))
  end

  weights = A\b
  return weights
end

function getWeights(domain::Ω, degree::Integer, points::AbstractVector{P}) where {Ω<:AbstractDomain,P<:AbstractVector}
  T = eltype(P)
  getWeights(T, domain, degree, points)
end