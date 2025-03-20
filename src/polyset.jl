# using PolyChaos: AbstractOrthoPoly

struct PolySet{Ω<:AbstractDomain}
  domain::Ω
  basis::Vector{Function}
end

struct BarycentricMonomials{dim, dim1} <: Function
  alpha::SVector{dim1,Int}

  function BarycentricMonomials(a::SVector{dim1,Int}) where {dim1}
    @assert length(a) > 0
    dim = length(a) - 1
    new{dim,dim1}(a)
  end
end


# λ -> λ[1]^alpha[1] * λ[2]^alpha[2] * ...
function (m::BarycentricMonomials{dim})(λ::AbstractVector) where {dim}
  @assert length(λ) == dim+1
  return prod(λ.^m.alpha)
end

# compute the integral of a barycentric monomial over a generic simplex of volume 1.
# For a specific simplex, this integral needs to be scaled with its volume.
function integrate(m::BarycentricMonomials{dim}, Ω::AbstractSimplex) where {dim}
  factorial(dim) * prod(factorial.(m.alpha)) // factorial(dim + sum(m.alpha))
end

function integrate(m::BarycentricMonomials{dim}, Ω::AbstractCube) where {dim}
  error("Not yet implemented")
end

function getTuples(length, total)
  if length == 1
    return NTuple{length,Int}[(total,)]
  end

  tuples = NTuple{length,Int}[]
  for i in 0:total
    for t in getTuples(length - 1, total - i)
      push!(tuples, (i,t...))
    end
  end
  return tuples
end

# Create a set of monomial of order <= degree
function PolySet(domain::Ω, degree::Integer) where {Ω<:AbstractDomain}
  dim = dimension(domain)

  basis = Function[]
  for deg in 0:degree
    for t in getTuples(dim+1, deg)
      push!(basis, BarycentricMonomials(SVector{dim+1,Int}(t...)))
    end
  end

  return PolySet(domain, basis)
end

struct LagrangeBasis{domain<:AbstractDomain}
  basis::Vector{Function}
  coeffs::Matrix{Rational}
  degree::Int
end

# using LinearAlgebra: I
# function LagrangeBasis(domain::Triangle, degree::Integer)
#   pointSet = EquidistantPointSet{Rational}(domain, degree)
#   polySet = PolySet(domain, degree)

#   @assert length(pointSet.points) == length(polySet.basis)

#   A = zeros(Rational, length(polyset.basis), length(pointSet.points))
#   for i in eachindex(polyset.basis)
#     f = polyset.basis[i]
#     for j in eachindex(pointSet.points)
#       A[i,j] = f(points[j])
#     end
#   end
#   B = I
#   coeffs = A\B

#   LagrangeBasis{Triangle}(polySet.basis, coeffs, degree)
# end

# function evaluateFunctions(basis::LagrangeBasis{Triangle}, λ::AbstractVector)
#   b = zeros(eltype(λ), length(basis.basis))
#   for i in eachindex(basis.basis)
#     f = basis.basis[i]
#     b[i] = f(λ)
#   end
#   basis.coeffs * b
# end

# function interpolate(basis::LagrangeBasis{Triangle}, f::Function)
#   pointSet = EquidistantPointSet{Rational}(Triangle(), basis.degree)
#   map(f, pointSet.points)
# end