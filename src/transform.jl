using StaticArrays: @SMatrix, @SVector, SMatrix, SVector

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


# transform the coordinates from the internal (barycentric) to the reference element domain
transformCoordinates(::D, X::AbstractVector{P}) where {D<:AbstractCube, P<:AbstractVector} = X
function transformCoordinates(::Triangle, X::AbstractVector{P}) where {P<:AbstractVector}
  let A = @SMatrix [-1 1 -1; -1 -1 1]
    map(λ -> A*λ, X)
  end
end
function transformCoordinates(::Tetrahedron, X::AbstractVector{P}) where {P<:AbstractVector}
  let A = @SMatrix [-1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1]
    map(λ -> A*λ, X)
  end
end
function transformCoordinates(::Prism, X::AbstractVector{P}) where {P<:AbstractVector}
  let A = SMatrix{3,4}(-1,-1,0, 1,-1,0, -1,1,0, 0,0,1),
      b = SVector{3}(0, 0, 0)
    map(x -> A*x + b, X)
  end
end
transformCoordinates(::Pyramid, X::AbstractVector{P}) where {P<:AbstractVector} = X


# transform the quadrature weights when changing the coordinates from barycentric to reference domain
transformWeights(::D, W::AbstractVector{<:Real}) where {D<:AbstractCube} = W
transformWeights(::Triangle, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Tetrahedron, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Prism, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Pyramid, W::AbstractVector{<:Real}) = W
