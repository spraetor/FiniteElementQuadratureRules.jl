using StaticArrays: SVector, MVector

"""
    LagrangePolySet{Ω<:AbstractDomain,R<:Real}

A `LagrangePolySet` represents a set of polynomials on a domain `Ω` in terms of a
set of basis polynomials, hereby given as Lagrange polynomials. As an `AbstractPolySet`
it provides also the integral values over the domain of all basis functions.
"""
struct LagrangePolySet{Ω<:AbstractDomain,R<:Real} <: AbstractPolySet
  domain::Ω
  basis::Vector{Function}
  integrals::Vector{R}
end


"""
    LagrangePolySet(::Type{T}, domain::AbstractDomain, degree::Integer)

Construct a `LagrangePolySet` on a given domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function LagrangePolySet(::Type{T}, domain::AbstractDomain, degree::Integer) where T
  basis = Function[]
  for k in 0:degree
    push!(basis, (x) -> _lagrange(domain, k, x))
  end
  integrals = zeros(T,length(basis)) # TODO: Not yet implemented.
  return LagrangePolySet(domain, basis, integrals)
end


"""
    LagrangePolySet(domain::AbstractDomain, degree::Integer)

Construct a `LagrangePolySet` on the given `domain` of given polynomial degree, with
`Float64` as data type used for the integral values.
"""
LagrangePolySet(domain::AbstractDomain, degree::Integer) = LagrangePolySet(Float64, domain, degree)


# evaluate all basis functions in a point x
function _lagrange(domain::AbstractSimplex, k::Integer, x::AbstractVector{T}) where {T<:Real}
  dim = dimension(domain)
  len = binomial(k+dim,dim)

  if k == 0
    SVector{1,T}(one(T))
  elseif k == 1
    barycentricCoordinates(domain, x)
  elseif k == 2
    λ = barycentricCoordinates(domain, x)
    if dim == 2
      SVector{len,T}(
        (λ.*(2 .* λ .- 1))...,             # vertex functions
        4λ[2]*λ[3], 4λ[3]*λ[1], 4λ[1]*λ[2] # edge functions
        )
    elseif dim == 3
      SVector{len,T}(
        (λ.*(2 .* λ .- 1))...,             # vertex functions
        4λ[1]*λ[2], 4λ[1]*λ[3], 4λ[1]*λ[4], 4λ[2]*λ[3], 4λ[2]*λ[4], 4λ[3]*λ[4] # edge functions
        )
    end
  elseif k == 3
    λ = barycentricCoordinates(domain, x)
    if dim == 2
      SVector{len,T}(
        ((3 .* λ .- 1).*(3 .* λ .- 2) .* λ ./ 2)...,    # vertex functions
        9(3λ[2]-1)*λ[2]*λ[3]/2, 9(3λ[3]-1)*λ[3]*λ[2]/2, # edge 1 functions
        9(3λ[3]-1)*λ[3]*λ[1]/2, 9(3λ[1]-1)*λ[1]*λ[3]/2, # edge 2 functions
        9(3λ[1]-1)*λ[1]*λ[2]/2, 9(3λ[2]-1)*λ[2]*λ[1]/2, # edge 3 functions
        27λ[1]*λ[2]*λ[3]                                # face functions
        )
    elseif dim == 3
      SVector{len,T}(
        ((3 .* λ .- 1).*(3 .* λ .- 2) .* λ ./ 2)...,    # vertex functions
        9(3λ[1]-1)*λ[1]*λ[2]/2, 9(3λ[2]-1)*λ[2]*λ[1]/2, # edge 1 functions
        9(3λ[1]-1)*λ[1]*λ[3]/2, 9(3λ[3]-1)*λ[3]*λ[1]/2, # edge 2 functions
        9(3λ[1]-1)*λ[1]*λ[4]/2, 9(3λ[4]-1)*λ[4]*λ[1]/2, # edge 3 functions
        9(3λ[2]-1)*λ[2]*λ[3]/2, 9(3λ[3]-1)*λ[3]*λ[2]/2, # edge 4 functions
        9(3λ[2]-1)*λ[2]*λ[4]/2, 9(3λ[4]-1)*λ[4]*λ[2]/2, # edge 5 functions
        9(3λ[3]-1)*λ[3]*λ[4]/2, 9(3λ[4]-1)*λ[4]*λ[3]/2, # edge 6 functions
        27λ[2]*λ[3]*λ[4], 27λ[3]*λ[4]*λ[1], 27λ[4]*λ[1]*λ[2], 27λ[1]*λ[2]*λ[3] # face functions
        )
    end
  else
    error("Not implemented")
  end
end


# check if the n'th bit is set in N
_bitisset(N::Int, n::Int) = (N >> (n-1)) & 1 != 0

# evaluate all basis functions in a point x
function _lagrange(domain::AbstractCube, k::Integer, x::AbstractVector{T}) where {T<:Real}
  dim = dimension(domain)
  len = (k+1)^dim

  if k == 0
    SVector{1,T}(one(T))
  elseif k == 1
    out = ones(T,len)
    for i in eachindex(out)
      for j in 1:dim
        out[i] *= _bitisset(i,j) ? (one(T)+x[j])/2 : (one(T)-x[j])/2
      end
    end
    SVector{len,T}(out)
  else
    error("Not implemented")
  end
end


# evaluate all basis functions in a point x
function _lagrange(::Prism, k::Integer, x::AbstractVector{T}) where {T<:Real}
  if k == 0
    SVector{1,T}(one(T))
  elseif k == 1
    λ = barycentricCoordinates(Triangle(), x[1:2])
    z = (one(T) + x[3]) / 2
    SVector{6,T}(
      λ[1]*(one(T)-z),
      λ[2]*(one(T)-z),
      λ[3]*(one(T)-z),
      λ[1]*z,
      λ[2]*z,
      λ[3]*z,
    )
  else
    error("Not implemented")
  end
end


# evaluate all basis functions in a point x
function _lagrange(::Pyramid, k::Integer, x::AbstractVector{T}) where {T<:Real}
  if k == 0
    SVector{1,T}(one(T))
  elseif k == 1
    ζ = (x[3] + one(T)) / 2
    d = one(T) - ζ

    # At the apex (z = 1), the four base shape functions vanish.
    if iszero(d)
      SVector{5,T}(zero(T), zero(T), zero(T), zero(T), one(T))
    else
      # Map from FEQ reference pyramid to the standard unit pyramid.
      ξ = (x[1] + d) / 2
      η = (x[2] + d) / 2

      SVector{5,T}(
        ((d-ξ)*(d-η))/d,
        (ξ*(d-η))/d,
        ((d-ξ)*η)/d,
        (ξ*η)/d,
        ζ
      )
    end
  else
    error("Not implemented")
  end
end
