# Transform a quadrature rule between reference elements

function transform(qr::QuadratureRule{T,D,Ω}, refIn::ReferenceElement{D,Ω}, refOut::ReferenceElement{D,Ω}) where {T,D,Ω}
  geo = Geometry(refIn, refOut.coordinates)
  volIn = volume(refIn)
  volOut = volume(refOut)

  QuadratureRule{T,D,Ω}(qr.degree,
    map(p -> map(geo,p), qr.points),
    map(w -> w*volOut/volIn, qr.weights),
    qr.properties, qr.accuracy, qr.bib)
end