using StaticArrays: SVector, SMatrix, MMatrix

"""
    AffineGeometry{Ω,mydim,cdom,T} <: AbstractGeometry

An `AffineGeometry` is a parametrization using an affine map of the form `x → A*x + b`.
As an `AbstractGeometry` it is itself a mapping.

Example:
```julia
geo = AffineGeometry(ReferenceElement(Triangle()), [x0,x1,x2])
@assert domaintype(geo) == Triangle()
x = geo(λ)    # given λ in coordinates of the reference element
```
"""
struct AffineGeometry{Ω<:AbstractDomain,mydim,cdim,T<:Real} <: AbstractGeometry
  refElement::ReferenceElement{Ω}
  origin::SVector{cdim,T}
  jacobian::SMatrix{cdim,mydim,T}
end


"""
    AffineGeometry(ref::ReferenceElement, coordVector::AbstractVector)

Given a reference element `ref` and a vector of corner coordinates `coordVector`,
construct an `AffineGeometry`.
"""
# constructor for the AffineGeometry
function AffineGeometry(ref::ReferenceElement{Ω}, coordVector::AbstractVector{C}) where {Ω<:AbstractDomain,cdim,T<:Real,C<:SVector{cdim,T}}
  mydim = dimension(Ω)
  @assert length(coordVector) == mydim+1

  jacobian = MMatrix{cdim,mydim,T}(undef)
  dref = MMatrix{mydim,mydim,T}(undef)
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

domaintype(::AffineGeometry{Ω,mydim,cdim,T}) where {Ω<:AbstractDomain,mydim,cdim,T} = Ω
domain(::AffineGeometry{Ω,mydim,cdim,T}) where {Ω<:AbstractDomain,mydim,cdim,T} = Ω()


"""
    (geo::AffineGeometry)(λ::AbstractVector)

Evaluate the affine geometry mapping `geo` in the reference element coordinate `λ`.
"""
function (geo::AffineGeometry)(λ::AbstractVector)
  geo.origin + geo.jacobian * λ
end


"""
    MultiLinearGeometry{T,cdim,Ω} <: AbstractGeometry

A `MultiLinearGeometry` is a geometry mapping defined in terms of linear Lagrange
basis function associated to the domain Ω. It can be constructed from a reference
element and a vector of corner coordinates.
"""
struct MultiLinearGeometry{Ω<:AbstractDomain,cdim,T<:Real} <: AbstractGeometry
  refElement::ReferenceElement{Ω}
  corners::Vector{SVector{cdim,T}}
end


"""
    MultiLinearGeometry(ref::ReferenceElement, coordVector::AbstractVector)

Construct a `MultiLinearGeometry` as an affine geometry if the domain is a simplex.
"""
function MultiLinearGeometry(ref::ReferenceElement{Ω}, coordVector::Vector{SVector{cdim,T}}) where {Ω<:AbstractSimplex,cdim,T<:Real}
  AffineGeometry(ref,coordVector)
end

domaintype(::MultiLinearGeometry{Ω,cdim,T}) where {Ω<:AbstractDomain,cdim,T} = Ω
domain(::MultiLinearGeometry{Ω,cdim,T}) where {Ω<:AbstractDomain,cdim,T} = Ω()

"""
    (geo::MultiLinearGeometry)(λ::AbstractVector)

Evaluate the multilinear geometry mapping `geo` in the reference element coordinates `λ`,
by linear combination of the associated Lagrange basis functions and corner vertices.
"""
function (geo::MultiLinearGeometry)(λ::AbstractVector)
  ϕs = _lagrange(domain(geo), 1, λ)

  @assert length(ϕs) == length(geo.corners)
  sum(ϕᵢ * cᵢ for (ϕᵢ,cᵢ) in zip(ϕs,geo.corners))
end
