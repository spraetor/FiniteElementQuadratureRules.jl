# Transform a quadrature rule between reference elements
function transform(qr::QuadratureRule{T,D,Ω}, refIn::ReferenceElement{D,Ω,P1}, refOut::ReferenceElement{D,Ω,P2}) where {T,D,Ω,P1,P2}
  geo = MultiLinearGeometry(refIn, refOut.coordinates)
  volIn = volume(refIn)
  volOut = volume(refOut)

  QuadratureRule{T,D,Ω}(qr.degree,
    map(geo, qr.points),
    map(w -> w*volOut/volIn, qr.weights),
    qr.properties, qr.accuracy, qr.bib)
end
