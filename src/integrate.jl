function integrate(f::Function, qr::QuadratureRule{Ω,T,P}) where {Ω<:AbstractDomain,T<:Real,P<:AbstractVector{T}}
  value::T = T(0)
  for i in 1:length(qr)
    value += f(qr.points[i]) * qr.weights[i]
  end
  return value
end