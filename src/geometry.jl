using StaticArrays: SVector, SMatrix, MMatrix

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

  jacobian = MMatrix{cdim,mydim,T}(zeros(T, cdim, mydim))
  dref = MMatrix{mydim,mydim,T}(zeros(T, mydim, mydim))
  for i in 1:mydim
    jacobian[:,i] .= coordVector[i+1] .- coordVector[1]
    dref[:,i] .= ref.coordinates[i+1] .- ref.coordinates[1]
  end

  # Reference coordinates are not necessarily the canonical unit simplex/cube coordinates.
  # Scale the Jacobian with the inverse reference edge matrix.
  if mydim > 0
    jacobian .= jacobian * inv(SMatrix(dref))
  end

  J = SMatrix(jacobian)
  x0 = mydim == 0 ? coordVector[1] : coordVector[1] - J * SVector{mydim,T}(ref.coordinates[1])

  AffineGeometry(ref, x0, J)
end

domaintype(::AffineGeometry{T,mydim,cdim,Ω}) where {T,mydim,cdim,Ω<:AbstractDomain} = Ω

# evaluate the affine geometry mapping
function (geo::AffineGeometry{T,mydim,cdim,Ω})(λ::AbstractVector{S}) where {T<:Real,mydim,cdim,Ω<:AbstractDomain,S<:Real}
  geo.origin + geo.jacobian * λ
end


"""
  MultiLinearGeometry{T,mydim,cdim,Ω}
"""
struct MultiLinearGeometry{T<:Real,mydim,cdim,Ω<:AbstractDomain} <: AbstractGeometry
  refElement::ReferenceElement{mydim,Ω}
  corners::Vector{SVector{cdim,T}}
end

# A multilinear geometry is an affine geometry if the reference element is a simplex
function MultiLinearGeometry(ref::ReferenceElement{mydim,Ω}, coordVector::Vector{SVector{cdim,T}}) where {mydim,Ω<:AbstractSimplex,cdim,T<:Real}
  AffineGeometry(ref,coordVector)
end

domaintype(::MultiLinearGeometry{T,mydim,cdim,Ω}) where {T,mydim,cdim,Ω<:AbstractDomain} = Ω


# evaluate the multilinear geometry mapping
function (geo::MultiLinearGeometry{T,mydim,cdim,Ω})(λ::AbstractVector{S}) where {T<:Real,mydim,cdim,Ω<:AbstractDomain,S<:Real}
  lb = LagrangeLocalBasis{mydim,Ω}(1)
  ϕs = lb(λ)

  @assert length(ϕs) == length(geo.corners)
  sum(ϕᵢ * cᵢ for (ϕᵢ,cᵢ) in zip(ϕs,geo.corners))
end
