using StaticArrays: @SMatrix, @SVector, SMatrix, SVector

# Transform a quadrature rule between reference elements
function transform(qr::QuadratureRule{Ω,T,P}, refIn::ReferenceElement{Ω,P1}, refOut::ReferenceElement{Ω,P2}) where {Ω,T,P,P1,P2}
  geo = MultiLinearGeometry(refIn, refOut.coordinates)
  volIn = volume(refIn)
  volOut = volume(refOut)

  QuadratureRule(qr.domain, qr.degree,
    map(geo, qr.points),
    map(w -> w*volOut/volIn, qr.weights))
end


# transform the coordinates from the internal (barycentric) to the reference element domain
transformCoordinates(::AbstractCube, X::AbstractVector{P}) where {P<:AbstractVector} = X
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

transformCoordinates(domain::AbstractDomain, X::NTuple{N,T}) where {N,T} = transformCoordinates(domain, SVector{N,T}(X))

# transform the quadrature weights when changing the coordinates from barycentric to reference domain
transformWeights(::AbstractCube, W::AbstractVector{<:Real}) = W
transformWeights(::Triangle, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Tetrahedron, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Prism, W::AbstractVector{<:Real}) = map(w -> 2*w, W)
transformWeights(::Pyramid, W::AbstractVector{<:Real}) = W
