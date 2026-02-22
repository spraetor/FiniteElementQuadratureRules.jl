using StaticArrays: SVector

"""
    JacobiPolySet{Ω<:AbstractDomain,R<:Real}

A `JacobiPolySet` represents a set of polynomials on a domain `Ω` in terms of a
set of basis polynomials, hereby given as Jacobi polynomials. As an `AbstractPolySet`
it provides also the integral values over the domain of all basis functions.
"""
struct JacobiPolySet{Ω<:AbstractDomain,R<:Real} <: AbstractPolySet
  domain::Ω
  basis::Vector{Function}
  integrals::Vector{R}
end


"""
    JacobiPolySet(::Type{T}, domain::Line, degree::Integer)

Construct a `JacobiPolySet` on the `Line` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Line, degree::Integer) where T
  basis = Function[]
  for i in 0:degree
    push!(basis, (x) -> _jacobi(i, x[1]))
  end

  integrals = zeros(T,length(basis))
  integrals[1] = 2
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(::Type{T}, domain::Triangle, degree::Integer)

Construct a `JacobiPolySet` on the `Triangle` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Triangle, degree::Integer) where T
  basis = Function[]
  for i in 0:degree
    for j in i:degree-i
      push!(basis, (x) -> let a = 2(1+x[1])/max(1-x[2], eps(T))-1,
                              b = x[2]
        sqrt(T(2))*_jacobi(i, a)*_jacobi(j,2i+1,0, b)*(1-b)^i
        end)
    end
  end

  integrals = zeros(T,length(basis))
  integrals[1] = 2*sqrt(T(2))
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(::Type{T}, domain::Quadrilateral, degree::Integer)

Construct a `JacobiPolySet` on the `Quadrilateral` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Quadrilateral, degree::Integer) where T
  basis = Function[]
  for i in 0:2:degree
    for j in i:2:degree-i
      push!(basis, (x) -> _jacobi(i, x[1])*_jacobi(j, x[2]))
    end
  end

  integrals = zeros(T,length(basis))
  integrals[1] = 4
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(::Type{T}, domain::Tetrahedron, degree::Integer)

Construct a `JacobiPolySet` on the `Tetrahedron` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Tetrahedron, degree::Integer) where T
  basis = Function[]
  for i in 0:degree
    for j in i:degree-i
      for k in j:degree-i-j
        push!(basis, (x) -> let a = -2(1+x[1])/(x[2]+x[3])-1,
                                b = 2(1+x[2])/max(1-x[3], eps(T))-1,
                                c = x[3]
          sqrt(T(8))*_jacobi(i, a)*_jacobi(j,2i+1,0, b)*_jacobi(k,2i+2j+2,0, c) * (1-b)^i * (1-c)^(i+j)
          end)
      end
    end
  end

  integrals = zeros(T,length(basis))
  integrals[1] = 4*sqrt((8))/T(3)
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(::Type{T}, domain::Prism, degree::Integer)

Construct a `JacobiPolySet` on the `Prism` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Prism, degree::Integer) where T
  basis = Function[]
  for i in 0:degree
    for j in i:degree-i
      for k in 0:2:degree-i-j
        push!(basis, (x) -> let a = 2(1+x[1])/max(1-x[2], eps(T))-1,
                                b = x[2],
                                c = x[3]
          sqrt(T(2))*_jacobi(i, a)*_jacobi(j,2i+1,0, b)*_jacobi(k, c)*(1-b)^i
          end)
      end
    end
  end

  integrals = zeros(T,length(basis))
  integrals[1] = 4*sqrt(T(2))
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(::Type{T}, domain::Pyramid, degree::Integer)

Construct a `JacobiPolySet` on the `Pyramid` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Pyramid, degree::Integer) where T
  basis = Function[]
  for i in 0:2:degree
    for j in i:2:degree-i
      for k in 0:degree-i-j
        push!(basis, (x) -> let a = 2x[1]/max(1-x[3], eps(T)),
                                b = 2x[2]/max(1-x[3], eps(T)),
                                c = x[3]
          2*_jacobi(i, a)*_jacobi(j, b)*_jacobi(k,2i+2j+2,0, c)*(1-c)^(i+j)
          end)
      end
    end
  end

  integrals = zeros(T,length(basis))
  integrals[1] = T(16)/T(3)
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(::Type{T}, domain::Hexahedron, degree::Integer)

Construct a `JacobiPolySet` on the `Hexahedron` domain of given polynomial degree, with
`T` the data type used for the integral values.
"""
function JacobiPolySet(::Type{T}, domain::Hexahedron, degree::Integer) where T
  basis = Function[]
  for i in 0:2:degree
    for j in i:2:degree-i
      for k in j:2:degree-i-j
        push!(basis, (x) -> _jacobi(i, x[1])*_jacobi(j, x[2])*_jacobi(k, x[3]))
      end
    end
  end

  integrals = zeros(T,length(basis))
  integrals[1] = 8
  return JacobiPolySet(domain, basis, integrals)
end


"""
    JacobiPolySet(domain::AbstractDomain, degree::Integer)

Construct a `JacobiPolySet` on the given `domain` of given polynomial degree, with
`Float64` as data type used for the integral values.
"""
JacobiPolySet(domain::AbstractDomain, degree::Integer) = JacobiPolySet(Float64, domain, degree)


# Evaluate all jacobi polynomials P_n^{a,b}(x) for n∈{0:N} and x∈X
function _jacobitable(N::Integer,a,b,X::AbstractVector)
  @assert N >= 0
  @assert a > -1 && b > -1

  P0 = ones(eltype(X),length(X))
  if N == 0
    return [P0]
  end

  P1 = collect((a .- b .+ (a + b + 2).*X)./2)
  if N == 1
    return [P0,P1]
  end

  P = zeros(eltype(X),(N+1,length(X)))
  P[1,:] .= P0
  P[2,:] .= P1
  for n in 1:N-1
    a1 = 2(n+1)*(n+a+b+1)*(2n+a+b)
    a2 = (2n+a+b+1)*(a^2-b^2)
    a3 = (2n+a+b)*(2n+a+b+1)*(2n+a+b+2)
    a4 = 2(n+a)*(n+b)*(2n+a+b+2)

    P0, P1 = P1, collect(((a2.+a3.*X).*P1 .- a4.*P0)./a1)
    P[n+2,:] .= P1
  end

  return P
end

# Evaluate the jacobi polynomial P_N^{a,b}(x) in all x∈X
function _jacobi(N::Integer,a,b,x)
  @assert N >= 0
  @assert a > -1 && b > -1

  P0 = one(x)
  if N == 0
    return P0
  end

  P1 = ((a - b) + (a + b + 2)*x)/2
  if N == 1
    return P1
  end

  P2 = zero(x)
  for n in 1:N-1
    a1 = 2(n+1)*(n+a+b+1)*(2n+a+b)
    a2 = (2n+a+b+1)*(a^2-b^2)
    a3 = (2n+a+b)*(2n+a+b+1)*(2n+a+b+2)
    a4 = 2(n+a)*(n+b)*(2n+a+b+2)

    P2 = ((a2 + a3*x)*P1 - a4*P0)/a1
    P0, P1 = P1, P2
  end

  return P2
end

# Evaluate the jacobi polynomial P_N^{a,b}(x) for all x∈X
function _jacobi(N::Integer,a,b,X::AbstractVector)
  map(x -> _jacobi(N,a,b,x), X)
end

# Evaluate the jacobi polynomial P_N^{0,0}(x) in x∈X
function _jacobi(N::Integer,x)
  @assert N >= 0

  P0 = one(x)
  if N == 0
    return P0
  end

  P1 = copy(x)
  if N == 1
    return P1
  end

  P2 = zero(x)
  for n in 1:N-1
    P2 = (((2n+1)*x)*P1 - n*P0)/(n+1)
    P0, P1 = P1, P2
  end

  return P2
end

# Evaluate the Jacobi polynomial P_N^{0,0}(x) in all x∈X
function _jacobi(N::Integer,X::AbstractVector)
  map(x -> _jacobi(N,x), X)
end

# Evaluate the derivative of the Jacobi polynomial, ∂ₓ P_N^{0,0}(x) in x∈X
function _djacobi(N::Integer,a,b,x)
  @assert N >= 0
  @assert a > -1 && b > -1

  if N == 0
    return zero(x)
  end

  if N == 1
    return (a+b+2)/2 * one(x)
  end

  P0 = one(x)
  P1 = (a - b + (a + b + 2)*x)/2
  P2 = zero(x)

  dP = zero(x)
  for n in 1:N-1
    a1 = 2(n+1)*(n+a+b+1)*(2n+a+b)
    a2 = (2n+a+b+1)*(a^2-b^2)
    a3 = (2n+a+b)*(2n+a+b+1)*(2n+a+b+2)
    a4 = 2(n+a)*(n+b)*(2n+a+b+2)

    P2 = ((a2+a3*x)*P1 - a4*P0)/a1
    P0, P1 = P1, P2

    b1 = (2n+a+b) * (1 - x^2)
    b2 = n((a - b) - (2n+a+b)*x)
    b3 = 2(n+a)*(n+b)

    dP = (b2*P1 + b3*P0)/b1
  end

  return dP
end

