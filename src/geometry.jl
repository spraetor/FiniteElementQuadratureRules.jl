using StaticArrays: SVector, SMatrix

abstract type AbstractGeometry end

struct AffineGeometry{T<:Real,mydim,cdim,Ω<:AbstractDomain} <: AbstractGeometry
  refElement::ReferenceElement{mydim,Ω}
  origin::SVector{cdim,T}
  jacobian::SMatrix{cdim,mydim,T}
end

function AffineGeometry(ref::ReferenceElement{mydim,Ω}, coordVector::AbstractVector{C}) where {mydim,Ω<:AbstractDomain,cdim,T<:Real,C<:SVector{cdim,T}}
  @assert length(coordVector) == mydim+1

  jacobian = MMatrix{cdim,mydim,T}()
  for i in 1:mydim
    jacobian[:,i] .= coordVector[i+1] .- coordVector[1]
  end

  AffineGeometry(ref,coordVector[1],SMatrix(jacobian))
end

import Base: map
function map(geo::AffineGeometry{T,mydim,cdim,Ω}, λ::AbstractVector{S}) where {mydim,Ω,cdim,T,S<:Real}
  geo.origin + geo.jacobian * λ
end