using StaticArrays: SVector

"""
    MonomialPolySet{Ω<:AbstractDomain,R<:Real}

A `MonomialPolySet` represents a set of polynomials on a domain `Ω` in terms of a
set of basis polynomials, hereby given as a monomial basis. As an `AbstractPolySet`
it provides also the integral values over the domain of all basis functions.
"""
struct MonomialPolySet{Ω<:AbstractDomain,R<:Real} <: AbstractPolySet
  domain::Ω
  basis::Vector{Function}
  integrals::Vector{R}
end

"""
    BarycentricMonomials <: Function

A representation of a multi-dimensional monomial basis function of the form
`λ[1]^α[1] * λ[2]^α[2] * ...` with `λ` the barycentric coordinate vector
and `α` the exponents associated to the components of the coordinates.
"""
struct BarycentricMonomials{dim, dim1} <: Function
  alpha::SVector{dim1,Int}

  function BarycentricMonomials(a::SVector{dim1,Int}) where {dim1}
    @assert length(a) > 0
    dim = length(a) - 1
    new{dim,dim1}(a)
  end
end


"""
    (m::BarycentricMonomials)(λ::AbstractVector)

Evaluate the monomial basis function `m` in a given barycentric coordinate `λ`, i.e.,
compute `λ -> λ[1]^α[1] * λ[2]^α[2] * ...`
"""
function (m::BarycentricMonomials{dim})(λ::AbstractVector) where {dim}
  @assert length(λ) == dim+1
  return prod(λ.^m.alpha)
end


"""
    integrate(m::BarycentricMonomials, Ω::AbstractSimplex)

Compute the integral of a barycentric monomial `m` over a generic simplex of volume 1.
For a specific simplex, this integral needs to be scaled with its volume.
"""
function integrate(m::BarycentricMonomials{dim}, Ω::AbstractSimplex) where {dim}
  factorial(dim) * prod(factorial.(m.alpha)) // factorial(dim + sum(m.alpha))
end

function integrate(m::BarycentricMonomials{dim}, Ω::AbstractDomain) where {dim}
  error("Not yet implemented")  # TODO
end


"""
    MonomialPolySet(::Type{T}, domain::AbstractSimplex, degree::Integer)

Construct a `MonomialPolySet` on a simplex domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function MonomialPolySet(::Type{T}, domain::AbstractSimplex, degree::Integer) where {T<:Real}
  dim = dimension(domain)

  basis = Function[]
  for deg in 0:degree
    for t in _allexponents(dim+1, deg)
      push!(basis, BarycentricMonomials(SVector{dim+1,Int}(t...)))
    end
  end

  integrals = map(b -> T(integrate(b,domain)), basis)
  return MonomialPolySet(domain, basis, integrals)
end

"""
    MonomialPolySet(domain::AbstractDomain, degree::Integer)

Construct a `MonomialPolySet` on the given `domain` of given polynomial degree, with
`Float64` as data type used for the integral values.
"""
MonomialPolySet(domain::AbstractDomain, degree::Integer) = MonomialPolySet(Float64,domain,degree)


# A utility to generate all possible exponents in a multidimensional polynomial
function _allexponents(length, total)
  if length == 1
    return NTuple{length,Int}[(total,)]
  end

  tuples = NTuple{length,Int}[]
  for i in 0:total
    for t in _allexponents(length - 1, total - i)
      push!(tuples, (i,t...))
    end
  end
  return tuples
end