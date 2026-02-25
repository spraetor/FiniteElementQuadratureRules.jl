using StaticArrays: SVector

struct SymmetryOrbit
  args::Union{Val{0},Val{1},Val{2},Val{3}}
  size::Int
  expand::Function
  compact::Function
end

import Base: length
length(so::SymmetryOrbit) = so.size

_args(::Val{0})::Int = 0
_args(::Val{1})::Int = 1
_args(::Val{2})::Int = 2
_args(::Val{3})::Int = 3

# Get the number of arguments
args(so::SymmetryOrbit)::Int = _args(so.args)

_expand(so::SymmetryOrbit, ::Val{0}, ::AbstractVector) = so.expand()
_expand(so::SymmetryOrbit, ::Val{1}, x::AbstractVector) = so.expand(x[1])
_expand(so::SymmetryOrbit, ::Val{2}, x::AbstractVector) = so.expand(x[1],x[2])
_expand(so::SymmetryOrbit, ::Val{3}, x::AbstractVector) = so.expand(x[1],x[2],x[3])

# Expand the symmetric orbit using the arguments
expand(so::SymmetryOrbit, x::AbstractVector) = _expand(so,so.args,x)

# Return the compact form of a coordinate
compact(so::SymmetryOrbit, x::AbstractVector) = so.compact(x)

function symmetryOrbits(Ω::AbstractDomain)
  symmetryOrbits(Float64, Ω)
end

function symmetryOrbits(::Type{T}, ::Point) where {T<:Real}
  P = SVector{1,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}(
      P(1),),
      (p)->NTuple{0,P}() )]
end

function symmetryOrbits(::Type{T}, ::Line) where {T<:Real}
  P = SVector{1,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}((
      P(0),)), (p)->SVector{1,P}() ),
    SymmetryOrbit(Val(1),2, (a::T) -> SVector{2,P}(
      P(a),
      P(-a)), (p)->(p[1],) )]
end

function symmetryOrbits(::Type{T}, ::Triangle) where {T<:Real}
  P = SVector{3,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}((            # S3
      P(1//3,1//3,1//3),)), (p)->NTuple{0,P}() ),
    SymmetryOrbit(Val(1),3, (a::T) -> SVector{3,P}(        # S21
      P(a,a,T(1)-2*a),
      P(a,T(1)-2*a,a),
      P(T(1)-2*a,a,a)), (p)->(p[1],) ),
    SymmetryOrbit(Val(2),6, (a::T, b::T) -> SVector{6,P}(  # S111
      P(a,b,T(1)-a-b),
      P(a,T(1)-a-b,b),
      P(T(1)-a-b,a,b),
      P(b,a,T(1)-a-b),
      P(b,T(1)-a-b,a),
      P(T(1)-a-b,b,a)), (p)->(p[1],p[2]) ),
    SymmetryOrbit(Val(2),2, (a::T, b::T) -> SVector{2,P}(  # T2=Mirror
      P(a,b,T(1)-a-b),
      P(T(1)-a-b,b,a)), (p)->(p[1],p[2]) ),
    SymmetryOrbit(Val(2),3, (a::T, b::T) -> SVector{3,P}(  # Ro3=Rotation=even permutations
      P(a,b,T(1)-a-b),
      P(b,T(1)-a-b,a),
      P(T(1)-a-b,a,b)), (p)->(p[1],p[2]) ), ]
end

function symmetryOrbits(::Type{T}, ::Quadrilateral) where {T<:Real}
  P = SVector{2,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}((
      P(0,0),)), (p)->NTuple{0,P}() ),
    SymmetryOrbit(Val(1),4, (a::T) -> SVector{4,P}(
      P(a,0),
      P(-a,0),
      P(0,a),
      P(0,-a)), (p)->!(p[1]≈0) ? (abs(p[1]),) : (abs(p[2]),) ),
    SymmetryOrbit(Val(1),4, (a::T) -> SVector{4,P}(
      P(a,a),
      P(-a,a),
      P(a,-a),
      P(-a,-a)), (p)->(abs(p[1]),) ),
    SymmetryOrbit(Val(2),8, (a::T, b::T) -> SVector{8,P}(
      P(a,b),
      P(b,a),
      P(-a,b),
      P(b,-a),
      P(a,-b),
      P(-b,a),
      P(-a,-b),
      P(-b,-a)), (p)->(abs(p[1]),abs(p[2])) )]
end

function symmetryOrbits(::Type{T}, ::Tetrahedron) where {T<:Real}
  P = SVector{4,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}((
      P(1//4,1//4,1//4,1//4),)), (p)->NTuple{0,P}() ),
    SymmetryOrbit(Val(1),4, (a::T) -> SVector{4,P}(
      P(a,a,a,T(1)-3*a),
      P(a,a,T(1)-3*a,a),
      P(a,T(1)-3*a,a,a),
      P(T(1)-3*a,a,a,a)), (p)->(p[1],) ),
    SymmetryOrbit(Val(1),6, (a::T) -> SVector{6,P}(
      P(a,a,1//2-a,1//2-a),
      P(a,1//2-a,a,1//2-a),
      P(1//2-a,a,a,1//2-a),
      P(a,1//2-a,1//2-a,a),
      P(1//2-a,a,1//2-a,a),
      P(1//2-a,1//2-a,a,a)), (p)->(p[1],) ),
    SymmetryOrbit(Val(2),12, (a::T, b::T) -> SVector{12,P}(
      P(a,a,b,T(1)-2*a-b),
      P(a,b,a,T(1)-2*a-b),
      P(b,a,a,T(1)-2*a-b),
      P(a,a,T(1)-2*a-b,b),
      P(a,T(1)-2*a-b,a,b),
      P(T(1)-2*a-b,a,a,b),
      P(a,b,T(1)-2*a-b,a),
      P(b,a,T(1)-2*a-b,a),
      P(a,T(1)-2*a-b,b,a),
      P(b,T(1)-2*a-b,a,a),
      P(T(1)-2*a-b,a,b,a),
      P(T(1)-2*a-b,b,a,a)), (p)->p[1]≈p[2] ? (p[1],p[3]) : (p[1],p[2]) ),
    SymmetryOrbit(Val(3),24, (a::T, b::T, c::T) -> SVector{24,P}(
      P(a,b,c,T(1)-a-b-c),
      P(a,c,b,T(1)-a-b-c),
      P(c,a,b,T(1)-a-b-c),
      P(c,b,a,T(1)-a-b-c),
      P(b,a,c,T(1)-a-b-c),
      P(b,c,a,T(1)-a-b-c),
      P(a,b,T(1)-a-b-c,c),
      P(a,c,T(1)-a-b-c,b),
      P(c,a,T(1)-a-b-c,b),
      P(c,b,T(1)-a-b-c,a),
      P(b,a,T(1)-a-b-c,c),
      P(b,c,T(1)-a-b-c,a),
      P(a,T(1)-a-b-c,b,c),
      P(a,T(1)-a-b-c,c,b),
      P(c,T(1)-a-b-c,a,b),
      P(c,T(1)-a-b-c,b,a),
      P(b,T(1)-a-b-c,a,c),
      P(b,T(1)-a-b-c,c,a),
      P(T(1)-a-b-c,a,b,c),
      P(T(1)-a-b-c,a,c,b),
      P(T(1)-a-b-c,c,a,b),
      P(T(1)-a-b-c,c,b,a),
      P(T(1)-a-b-c,b,a,c),
      P(T(1)-a-b-c,b,c,a)), (p)->(p[1],p[2],p[3]) )]
end

function symmetryOrbits(::Type{T}, ::Hexahedron) where {T<:Real}
  P = SVector{3,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}((
      P(0,0,0),)), (p)->NTuple{0,P}() ),
    SymmetryOrbit(Val(1),6, (a::T) -> SVector{6,P}(
      P(a,0,0),
      P(-a,0,0),
      P(0,a,0),
      P(0,-a,0),
      P(0,0,a),
      P(0,0,-a)), (p)->!(p[1]≈0) ? (abs(p[1]),) : !(p[2]≈0) ? (abs(p[2]),) : (abs(p[3]),) ),
    SymmetryOrbit(Val(1),8, (a::T) -> SVector{8,P}(
      P(a,a,a),
      P(-a,a,a),
      P(a,-a,a),
      P(a,a,-a),
      P(-a,-a,a),
      P(-a,a,-a),
      P(a,-a,-a),
      P(-a,-a,-a)), (p)->(abs(p[1]),) ),
    SymmetryOrbit(Val(1),12, (a::T) -> SVector{12,P}(
      P(a,a,0),
      P(a,0,a),
      P(0,a,a),
      P(a,-a,0),
      P(a,0,-a),
      P(0,a,-a),
      P(-a,a,0),
      P(-a,0,a),
      P(0,-a,a),
      P(-a,-a,0),
      P(-a,0,-a),
      P(0,-a,-a)), (p)->!(p[1]≈0) ? (abs(p[1]),) : !(p[2]≈0) ? (abs(p[2]),) : (abs(p[3]),) ),
    SymmetryOrbit(Val(2),24, (a::T, b::T) -> SVector{24,P}(
      P(a,b,0),
      P(a,0,b),
      P(b,a,0),
      P(b,0,a),
      P(0,a,b),
      P(0,b,a),
      P(-a,b,0),
      P(-a,0,b),
      P(b,-a,0),
      P(b,0,-a),
      P(0,-a,b),
      P(0,b,-a),
      P(a,-b,0),
      P(a,0,-b),
      P(-b,a,0),
      P(-b,0,a),
      P(0,a,-b),
      P(0,-b,a),
      P(-a,-b,0),
      P(-a,0,-b),
      P(-b,-a,0),
      P(-b,0,-a),
      P(0,-a,-b),
      P(0,-b,-a)), (p)->p[1]≈0 ? (abs(p[2]),abs(p[3])) : p[2]≈0 ? (abs(p[1]),abs(p[2])) : (abs(p[1]),abs(p[2])) ),
    SymmetryOrbit(Val(2),24, (a::T, b::T) -> SVector{24,P}(
      P(a,a,b),
      P(a,b,a),
      P(b,a,a),
      P(-a,a,b),
      P(-a,b,a),
      P(a,-a,b),
      P(b,-a,a),
      P(a,b,-a),
      P(b,a,-a),
      P(-a,-a,b),
      P(-a,b,-a),
      P(b,-a,-a),
      P(-a,a,-b),
      P(-a,-b,a),
      P(a,-a,-b),
      P(-b,-a,a),
      P(a,-b,-a),
      P(-b,a,-a),
      P(a,a,-b),
      P(a,-b,a),
      P(-b,a,a),
      P(-a,-a,-b),
      P(-a,-b,-a),
      P(-b,-a,-a)), (p)->p[1]≈p[2] ? (abs(p[1]),abs(p[3])) : p[1]≈p[3] ?  (abs(p[1]),abs(p[2])) :  (abs(p[2]),abs(p[3])) ),
    SymmetryOrbit(Val(3),48, (a::T, b::T, c::T) -> SVector{48,P}(
      P(a,b,c),
      P(a,c,b),
      P(b,a,c),
      P(b,c,a),
      P(c,a,b),
      P(c,b,a),
      P(-a,b,c),
      P(-a,c,b),
      P(b,-a,c),
      P(b,c,-a),
      P(c,-a,b),
      P(c,b,-a),
      P(a,-b,c),
      P(a,c,-b),
      P(-b,a,c),
      P(-b,c,a),
      P(c,a,-b),
      P(c,-b,a),
      P(a,b,-c),
      P(a,-c,b),
      P(b,a,-c),
      P(b,-c,a),
      P(-c,a,b),
      P(-c,b,a),
      P(-a,-b,c),
      P(-a,c,-b),
      P(-b,-a,c),
      P(-b,c,-a),
      P(c,-a,-b),
      P(c,-b,-a),
      P(-a,b,-c),
      P(-a,-c,b),
      P(b,-a,-c),
      P(b,-c,-a),
      P(-c,-a,b),
      P(-c,b,-a),
      P(a,-b,-c),
      P(a,-c,-b),
      P(-b,a,-c),
      P(-b,-c,a),
      P(-c,a,-b),
      P(-c,-b,a),
      P(-a,-b,-c),
      P(-a,-c,-b),
      P(-b,-a,-c),
      P(-b,-c,-a),
      P(-c,-a,-b),
      P(-c,-b,-a)), (p)->(abs(p[1]),abs(p[2]),abs(p[3])) )]
end

function symmetryOrbits(::Type{T}, ::Prism) where {T<:Real}
  P = SVector{4,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(0),1, () -> SVector{1,P}((
      P(1//3,1//3,1//3,0),)), (p)->NTuple{0,P}() ),
    SymmetryOrbit(Val(1),2, (c::T) -> SVector{2,P}(
      P(1//3,1//3,1//3,c),
      P(1//3,1//3,1//3,-c)), (p)->(abs(p[4]),) ),
    SymmetryOrbit(Val(1),3, (a::T) -> SVector{3,P}(
      P(a,a,T(1)-2*a,0),
      P(a,T(1)-2*a,a,0),
      P(T(1)-2*a,a,a,0)), (p)->(p[1],) ),
    SymmetryOrbit(Val(2),6, (a::T, c::T) -> SVector{6,P}(
      P(a,a,T(1)-2*a,c),
      P(a,T(1)-2*a,a,c),
      P(T(1)-2*a,a,a,c),
      P(a,a,T(1)-2*a,-c),
      P(a,T(1)-2*a,a,-c),
      P(T(1)-2*a,a,a,-c)), (p)->(p[1],abs(p[4])) ),
    SymmetryOrbit(Val(2),6, (a::T, b::T) -> SVector{6,P}(
      P(a,b,T(1)-a-b,0),
      P(a,T(1)-a-b,b,0),
      P(T(1)-a-b,a,b,0),
      P(b,a,T(1)-a-b,0),
      P(b,T(1)-a-b,a,0),
      P(T(1)-a-b,b,a,0)), (p)->(p[1],p[2]) ),
    SymmetryOrbit(Val(3),12, (a::T, b::T, c::T) -> SVector{12,P}(
      P(a,b,T(1)-a-b,c),
      P(a,T(1)-a-b,b,c),
      P(T(1)-a-b,a,b,c),
      P(b,a,T(1)-a-b,c),
      P(b,T(1)-a-b,a,c),
      P(T(1)-a-b,b,a,c),
      P(a,b,T(1)-a-b,-c),
      P(a,T(1)-a-b,b,-c),
      P(T(1)-a-b,a,b,-c),
      P(b,a,T(1)-a-b,-c),
      P(b,T(1)-a-b,a,-c),
      P(T(1)-a-b,b,a,-c)), (p)->(p[1],p[2],abs(p[4])) )]
end

function symmetryOrbits(::Type{T}, ::Pyramid) where {T<:Real}
  P = SVector{3,T}
  SymmetryOrbit[
    SymmetryOrbit(Val(1),1, (c::T) -> SVector{1,P}((
      P(0,0,c),)), (p)->(p[3],) ),
    SymmetryOrbit(Val(2),4, (a::T, c::T) -> SVector{4,P}(
      P(a,0,c),
      P(-a,0,c),
      P(0,a,c),
      P(0,-a,c)), (p)->p[1] != 0 ? (abs(p[1]),p[3]) : (abs(p[2]),p[3]) ),
    SymmetryOrbit(Val(2),4, (a::T, c::T) -> SVector{4,P}(
      P(a,a,c),
      P(-a,a,c),
      P(a,-a,c),
      P(-a,-a,c)), (p)->(abs(p[1]),p[3]) ),
    SymmetryOrbit(Val(3),8, (a::T, b::T, c::T) -> SVector{8,P}(
      P(a,b,c),
      P(b,a,c),
      P(-a,b,c),
      P(b,-a,c),
      P(a,-b,c),
      P(-b,a,c),
      P(-a,-b,c),
      P(-b,-a,c)), (p)->(abs(p[1]),abs(p[2]),p[3]) )]
end
