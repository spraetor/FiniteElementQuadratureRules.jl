using StaticArrays: SVector, SMatrix

"""
  AffineGeometry{T,mydim,cdom,Ω}
"""
struct AffineGeometry{T<:Real,mydim,cdim,Ω<:AbstractDomain} <: AbstractGeometry
  refElement::ReferenceElement{mydim,Ω}
  origin::SVector{cdim,T}
  jacobian::SMatrix{cdim,mydim,T}
end

# constructor for the AffineGeometry
function AffineGeometry(ref::ReferenceElement{mydim,Ω}, coordVector::AbstractVector{C}) where {mydim,Ω<:AbstractDomain,cdim,T<:Real,C<:SVector{cdim,T}}
  @assert length(coordVector) == mydim+1

  jacobian = MMatrix{cdim,mydim,T}()
  for i in 1:mydim
    jacobian[:,i] .= coordVector[i+1] .- coordVector[1]
  end

  AffineGeometry(ref,coordVector[1],SMatrix(jacobian))
end

# evaluate the affine geometry mapping
function (geo::AffineGeometry{T,mydim,cdim,Ω})(λ::AbstractVector{S}) where {T<:Real,mydim,cdim,Ω<:AbstractDomain,S<:Real}
  geo.origin + geo.jacobian * λ
end


"""
  MultiLinearGeometry{T,mydim,cdom,Ω}
"""
struct MultiLinearGeometry{T<:Real,mydim,cdim,Ω<:AbstractDomain} <: AbstractGeometry
  refElement::ReferenceElement{mydim,Ω}
  corners::Vector{SVector{cdim,T}}
end

# A multilinear geometry is an affine geometry if the reference element is a simplex
function MultiLinearGeometry(ref::ReferenceElement{mydim,Ω}, coordVector::AbstractVector{C}) where {mydim,Ω<:AbstractSimplex,cdim,T<:Real,C<:SVector{cdim,T}}
  AffineGeometry(ref,coordVector)
end

# evaluate the multilinear geometry mapping
function (geo::MultiLinearGeometry{T,mydim,cdim,Ω})(λ::AbstractVector{S}) where {T<:Real,mydim,cdim,Ω<:AbstractDomain,S<:Real}
  lb = Lagrange{mydim,Ω}(1)
  ϕs = lb(λ)

  @assert length(ϕs) == length(geo.corners)
  sum(ϕᵢ * cᵢ for (ϕᵢ,cᵢ) in zip(ϕs,geo.corners))
end