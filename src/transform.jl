# Transform a quadrature rule between reference elements
function transform(qr::QuadratureRule{Ω,T,P}, refOut::ReferenceElement) where {Ω,T,P}
  @assert domain(qr) == domain(refOut)
  if qr.ref == refOut
    return qr
  end

  Ref = typeof(qr.ref)
  refOutT = Ref(refOut)
  geo = MultiLinearGeometry(qr.ref, coordinates(refOutT))
  volIn = volume(qr.ref)
  volOut = volume(refOutT)

  QuadratureRule(refOutT, qr.degree,
    map(geo, qr.points),
    map(w -> w*volOut/volIn, qr.weights),
    qr.properties)
end

