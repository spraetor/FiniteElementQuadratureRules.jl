using StaticArrays: SVector

"""
  Lagrange{dim,Ω}
"""
struct Lagrange{dim,Ω<:AbstractDomain} <: AbstractLocalBasis
  k::Int
end

import Base: length
length(lb::Lagrange{dim,<:AbstractSimplex}) where dim = binomial(lb.k+dim,dim)
length(lb::Lagrange{dim,<:AbstractCube}) where dim = (lb.k+1)^dim

# Return the polynomial order of the basis
order(lb::Lagrange) = lb.k

# Return the domain type
domain(::Lagrange{dim,Ω}) where {dim,Ω<:AbstractSimplex} = Ω()


# evaluate all basis functions in a point x
function (lb::Lagrange{dim,<:AbstractSimplex})(x::AbstractVector{T}) where {dim,T<:Real}
  if order(lb) == 0
    SVector{1,T}([1.0])
  elseif order(lb) == 1
    SVector{dim+1,T}([1.0 - sum(x), x])
  else
    error("Not implemented")
  end
end


# check if the n'th bit is set in N
bit_is_set(N::Int, n::Int) = (N >> (n-1)) & 1 != 0

# evaluate all basis functions in a point x
function (lb::Lagrange{dim,<:AbstractCube})(x::AbstractVector{T}) where {dim,T<:Real}
  if order(lb) == 0
    SVector{1,T}([1.0])
  elseif order(lb) == 1
    out = MVector{2^dim,T}(ones(T,length(lb)))
    for i in eachindex(out)
      for j in 1:dim
        out[i] *= bit_is_set(i,j) ? x[j] : 1-x[j];
      end
    end
    SVector(out)
  else
    error("Not implemented")
  end
end


# evaluate all basis functions in a point x
function (lb::Lagrange{3,Prism})(x::AbstractVector{T}) where {dim,T<:Real}
  if order(lb) == 0
    SVector{1,T}([1.0])
  elseif order(lb) == 1
    SVector{6,T}([
      (1.0-x[1]-x[2])*(1.0-x[3]),
      x[1]*(1-x[3]),
      x[2]*(1-x[3]),
      x[3]*(1.0-x[1]-x[2]),
      x[1]*x[3],
      x[2]*x[3]
    ])
  else
    error("Not implemented")
  end
end


# evaluate all basis functions in a point x
function (lb::Lagrange{3,Pyramid})(x::AbstractVector{T}) where {dim,T<:Real}
  if order(lb) == 0
    SVector{1,T}([1.0])
  elseif order(lb) == 1
    if x[1] > x[2]
      SVector{5,T}([
        (1-x[1])*(1-x[2])-x[3]*(1-x[2]),
        x[1]*(1-x[2])-x[3]*x[2],
        (1-x[1])*x[2]-x[3]*x[2],
        x[1]*x[2]+x[3]*x[2],
        x[3]
      ])
    else
      SVector{5,T}([
        (1-x[1])*(1-x[2])-x[3]*(1-x[1]),
        x[1]*(1-x[2])-x[3]*x[1],
        (1-x[1])*x[2]-x[3]*x[1],
        x[1]*x[2]+x[3]*x[1],
        x[3]
      ])
    end
  else
    error("Not implemented")
  end
end