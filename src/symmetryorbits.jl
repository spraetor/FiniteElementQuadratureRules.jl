using StaticArrays: SVector

struct SymmetryOrbit
  args::Int
  size::Int
  orbit::Function
end

import Base: length
length(so::SymmetryOrbit) = so.size

function symmetryOrbits(Ω::AbstractDomain)
  symmetryOrbits(Float64, Ω)
end

function symmetryOrbits(::Type{T}, ::Point) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (
      SVector{1,T}(1),) )]
end

function symmetryOrbits(::Type{T}, ::Line) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (
      SVector{1,T}(0),) ),
    SymmetryOrbit(1,2, (a::T) -> (
      SVector{1,T}(a),
      SVector{1,T}(-a)) )]
end

function symmetryOrbits(::Type{T}, ::Triangle) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (            # S3
      SVector{3,T}(1//3,1//3,1//3),) ),
    SymmetryOrbit(1,3, (a::T) -> (        # S21
      SVector{3,T}(a,a,T(1)-2*a),
      SVector{3,T}(a,T(1)-2*a,a),
      SVector{3,T}(T(1)-2*a,a,a)) ),
    SymmetryOrbit(2,6, (a::T, b::T) -> (  # S111
      SVector{3,T}(a,b,T(1)-a-b),
      SVector{3,T}(a,T(1)-a-b,b),
      SVector{3,T}(T(1)-a-b,a,b),
      SVector{3,T}(b,a,T(1)-a-b),
      SVector{3,T}(b,T(1)-a-b,a),
      SVector{3,T}(T(1)-a-b,b,a)) )]
end

function symmetryOrbits(::Type{T}, ::Quadrilateral) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (
      SVector{2,T}(0,0),) ),
    SymmetryOrbit(1,4, (a::T) -> (
      SVector{2,T}(a,0),
      SVector{2,T}(-a,0),
      SVector{2,T}(0,a),
      SVector{2,T}(0,-a)) ),
    SymmetryOrbit(1,4, (a::T) -> (
      SVector{2,T}(a,a),
      SVector{2,T}(-a,a),
      SVector{2,T}(a,-a),
      SVector{2,T}(-a,-a)) ),
    SymmetryOrbit(2,8, (a::T, b::T) -> (
      SVector{2,T}(a,b),
      SVector{2,T}(b,a),
      SVector{2,T}(-a,b),
      SVector{2,T}(b,-a),
      SVector{2,T}(a,-b),
      SVector{2,T}(-b,a),
      SVector{2,T}(-a,-b),
      SVector{2,T}(-b,-a)) )]
end

function symmetryOrbits(::Type{T}, ::Tetrahedron) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (
      SVector{4,T}(1//4,1//4,1//4,1//4),) ),
    SymmetryOrbit(1,4, (a::T) -> (
      SVector{4,T}(a,a,a,T(1)-3*a),
      SVector{4,T}(a,a,T(1)-3*a,a),
      SVector{4,T}(a,T(1)-3*a,a,a),
      SVector{4,T}(T(1)-3*a,a,a,a)) ),
    SymmetryOrbit(1,6, (a::T) -> (
      SVector{4,T}(a,a,1//2-a,1//2-a),
      SVector{4,T}(a,1//2-a,a,1//2-a),
      SVector{4,T}(1//2-a,a,a,1//2-a),
      SVector{4,T}(a,1//2-a,1//2-a,a),
      SVector{4,T}(1//2-a,a,1//2-a,a),
      SVector{4,T}(1//2-a,1//2-a,a,a)) ),
    SymmetryOrbit(2,12, (a::T, b::T) -> (
      SVector{4,T}(a,a,b,T(1)-2*a-b),
      SVector{4,T}(a,b,a,T(1)-2*a-b),
      SVector{4,T}(b,a,a,T(1)-2*a-b),
      SVector{4,T}(a,a,T(1)-2*a-b,b),
      SVector{4,T}(a,T(1)-2*a-b,a,b),
      SVector{4,T}(T(1)-2*a-b,a,a,b),
      SVector{4,T}(a,b,T(1)-2*a-b,a),
      SVector{4,T}(b,a,T(1)-2*a-b,a),
      SVector{4,T}(a,T(1)-2*a-b,b,a),
      SVector{4,T}(b,T(1)-2*a-b,a,a),
      SVector{4,T}(T(1)-2*a-b,a,b,a),
      SVector{4,T}(T(1)-2*a-b,b,a,a)) ),
    SymmetryOrbit(3,24, (a::T, b::T, c::T) -> (
      SVector{4,T}(a,b,c,T(1)-a-b-c),
      SVector{4,T}(a,c,b,T(1)-a-b-c),
      SVector{4,T}(c,a,b,T(1)-a-b-c),
      SVector{4,T}(c,b,a,T(1)-a-b-c),
      SVector{4,T}(b,a,c,T(1)-a-b-c),
      SVector{4,T}(b,c,a,T(1)-a-b-c),
      SVector{4,T}(a,b,T(1)-a-b-c,c),
      SVector{4,T}(a,c,T(1)-a-b-c,b),
      SVector{4,T}(c,a,T(1)-a-b-c,b),
      SVector{4,T}(c,b,T(1)-a-b-c,a),
      SVector{4,T}(b,a,T(1)-a-b-c,c),
      SVector{4,T}(b,c,T(1)-a-b-c,a),
      SVector{4,T}(a,T(1)-a-b-c,b,c),
      SVector{4,T}(a,T(1)-a-b-c,c,b),
      SVector{4,T}(c,T(1)-a-b-c,a,b),
      SVector{4,T}(c,T(1)-a-b-c,b,a),
      SVector{4,T}(b,T(1)-a-b-c,a,c),
      SVector{4,T}(b,T(1)-a-b-c,c,a),
      SVector{4,T}(T(1)-a-b-c,a,b,c),
      SVector{4,T}(T(1)-a-b-c,a,c,b),
      SVector{4,T}(T(1)-a-b-c,c,a,b),
      SVector{4,T}(T(1)-a-b-c,c,b,a),
      SVector{4,T}(T(1)-a-b-c,b,a,c),
      SVector{4,T}(T(1)-a-b-c,b,c,a)) )]
end

function symmetryOrbits(::Type{T}, ::Hexahedron) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (
      SVector{3,T}(0,0,0),) ),
    SymmetryOrbit(1,6, (a::T) -> (
      SVector{3,T}(a,0,0),
      SVector{3,T}(-a,0,0),
      SVector{3,T}(0,a,0),
      SVector{3,T}(0,-a,0),
      SVector{3,T}(0,0,a),
      SVector{3,T}(0,0,-a)) ),
    SymmetryOrbit(1,8, (a::T) -> (
      SVector{3,T}(a,a,a),
      SVector{3,T}(-a,a,a),
      SVector{3,T}(a,-a,a),
      SVector{3,T}(a,a,-a),
      SVector{3,T}(-a,-a,a),
      SVector{3,T}(-a,a,-a),
      SVector{3,T}(a,-a,-a),
      SVector{3,T}(-a,-a,-a)) ),
    SymmetryOrbit(1,12, (a::T) -> (
      SVector{3,T}(a,a,0),
      SVector{3,T}(a,0,a),
      SVector{3,T}(0,a,a),
      SVector{3,T}(a,-a,0),
      SVector{3,T}(a,0,-a),
      SVector{3,T}(0,a,-a),
      SVector{3,T}(-a,a,0),
      SVector{3,T}(-a,0,a),
      SVector{3,T}(0,-a,a),
      SVector{3,T}(-a,-a,0),
      SVector{3,T}(-a,0,-a),
      SVector{3,T}(0,-a,-a)) ),
    SymmetryOrbit(2,24, (a::T, b::T) -> (
      SVector{3,T}(a,b,0),
      SVector{3,T}(a,0,b),
      SVector{3,T}(b,a,0),
      SVector{3,T}(b,0,a),
      SVector{3,T}(0,a,b),
      SVector{3,T}(0,b,a),
      SVector{3,T}(-a,b,0),
      SVector{3,T}(-a,0,b),
      SVector{3,T}(b,-a,0),
      SVector{3,T}(b,0,-a),
      SVector{3,T}(0,-a,b),
      SVector{3,T}(0,b,-a),
      SVector{3,T}(a,-b,0),
      SVector{3,T}(a,0,-b),
      SVector{3,T}(-b,a,0),
      SVector{3,T}(-b,0,a),
      SVector{3,T}(0,a,-b),
      SVector{3,T}(0,-b,a),
      SVector{3,T}(-a,-b,0),
      SVector{3,T}(-a,0,-b),
      SVector{3,T}(-b,-a,0),
      SVector{3,T}(-b,0,-a),
      SVector{3,T}(0,-a,-b),
      SVector{3,T}(0,-b,-a)) ),
    SymmetryOrbit(2,24, (a::T, b::T) -> (
      SVector{3,T}(a,a,b),
      SVector{3,T}(a,b,a),
      SVector{3,T}(b,a,a),
      SVector{3,T}(-a,a,b),
      SVector{3,T}(-a,b,a),
      SVector{3,T}(a,-a,b),
      SVector{3,T}(b,-a,a),
      SVector{3,T}(a,b,-a),
      SVector{3,T}(b,a,-a),
      SVector{3,T}(-a,-a,b),
      SVector{3,T}(-a,b,-a),
      SVector{3,T}(b,-a,-a),
      SVector{3,T}(-a,a,-b),
      SVector{3,T}(-a,-b,a),
      SVector{3,T}(a,-a,-b),
      SVector{3,T}(-b,-a,a),
      SVector{3,T}(a,-b,-a),
      SVector{3,T}(-b,a,-a),
      SVector{3,T}(a,a,-b),
      SVector{3,T}(a,-b,a),
      SVector{3,T}(-b,a,a),
      SVector{3,T}(-a,-a,-b),
      SVector{3,T}(-a,-b,-a),
      SVector{3,T}(-b,-a,-a)) ),
    SymmetryOrbit(3,48, (a::T, b::T, c::T) -> (
      SVector{3,T}(a,b,c),
      SVector{3,T}(a,c,b),
      SVector{3,T}(b,a,c),
      SVector{3,T}(b,c,a),
      SVector{3,T}(c,a,b),
      SVector{3,T}(c,b,a),
      SVector{3,T}(-a,b,c),
      SVector{3,T}(-a,c,b),
      SVector{3,T}(b,-a,c),
      SVector{3,T}(b,c,-a),
      SVector{3,T}(c,-a,b),
      SVector{3,T}(c,b,-a),
      SVector{3,T}(a,-b,c),
      SVector{3,T}(a,c,-b),
      SVector{3,T}(-b,a,c),
      SVector{3,T}(-b,c,a),
      SVector{3,T}(c,a,-b),
      SVector{3,T}(c,-b,a),
      SVector{3,T}(a,b,-c),
      SVector{3,T}(a,-c,b),
      SVector{3,T}(b,a,-c),
      SVector{3,T}(b,-c,a),
      SVector{3,T}(-c,a,b),
      SVector{3,T}(-c,b,a),
      SVector{3,T}(-a,-b,c),
      SVector{3,T}(-a,c,-b),
      SVector{3,T}(-b,-a,c),
      SVector{3,T}(-b,c,-a),
      SVector{3,T}(c,-a,-b),
      SVector{3,T}(c,-b,-a),
      SVector{3,T}(-a,b,-c),
      SVector{3,T}(-a,-c,b),
      SVector{3,T}(b,-a,-c),
      SVector{3,T}(b,-c,-a),
      SVector{3,T}(-c,-a,b),
      SVector{3,T}(-c,b,-a),
      SVector{3,T}(a,-b,-c),
      SVector{3,T}(a,-c,-b),
      SVector{3,T}(-b,a,-c),
      SVector{3,T}(-b,-c,a),
      SVector{3,T}(-c,a,-b),
      SVector{3,T}(-c,-b,a),
      SVector{3,T}(-a,-b,-c),
      SVector{3,T}(-a,-c,-b),
      SVector{3,T}(-b,-a,-c),
      SVector{3,T}(-b,-c,-a),
      SVector{3,T}(-c,-a,-b),
      SVector{3,T}(-c,-b,-a)) )]
end

function symmetryOrbits(::Type{T}, ::Prism) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(0,1, () -> (
      SVector{4,T}(1//3,1//3,1//3,0),) ),
    SymmetryOrbit(1,2, (c::T) -> (
      SVector{4,T}(1//3,1//3,1//3,c),
      SVector{4,T}(1//3,1//3,1//3,-c)) ),
    SymmetryOrbit(1,3, (a::T) -> (
      SVector{4,T}(a,a,T(1)-2*a,0),
      SVector{4,T}(a,T(1)-2*a,a,0),
      SVector{4,T}(T(1)-2*a,a,a,0)) ),
    SymmetryOrbit(2,6, (a::T, c::T) -> (
      SVector{4,T}(a,a,T(1)-2*a,c),
      SVector{4,T}(a,T(1)-2*a,a,c),
      SVector{4,T}(T(1)-2*a,a,a,c),
      SVector{4,T}(a,a,T(1)-2*a,-c),
      SVector{4,T}(a,T(1)-2*a,a,-c),
      SVector{4,T}(T(1)-2*a,a,a,-c)) ),
    SymmetryOrbit(2,6, (a::T, b::T) -> (
      SVector{4,T}(a,b,T(1)-a-b,0),
      SVector{4,T}(a,T(1)-a-b,b,0),
      SVector{4,T}(T(1)-a-b,a,b,0),
      SVector{4,T}(b,a,T(1)-a-b,0),
      SVector{4,T}(b,T(1)-a-b,a,0),
      SVector{4,T}(T(1)-a-b,b,a,0)) ),
    SymmetryOrbit(3,12, (a::T, b::T, c::T) -> (
      SVector{4,T}(a,b,T(1)-a-b,c),
      SVector{4,T}(a,T(1)-a-b,b,c),
      SVector{4,T}(T(1)-a-b,a,b,c),
      SVector{4,T}(b,a,T(1)-a-b,c),
      SVector{4,T}(b,T(1)-a-b,a,c),
      SVector{4,T}(T(1)-a-b,b,a,c),
      SVector{4,T}(a,b,T(1)-a-b,-c),
      SVector{4,T}(a,T(1)-a-b,b,-c),
      SVector{4,T}(T(1)-a-b,a,b,-c),
      SVector{4,T}(b,a,T(1)-a-b,-c),
      SVector{4,T}(b,T(1)-a-b,a,-c),
      SVector{4,T}(T(1)-a-b,b,a,-c)) )]
end

function symmetryOrbits(::Type{T}, ::Pyramid) where {T<:Real}
  SymmetryOrbit[
    SymmetryOrbit(1,1, (c::T) -> (
      SVector{3,T}(0,0,c),) ),
    SymmetryOrbit(2,4, (a::T, c::T) -> (
      SVector{3,T}(a,0,c),
      SVector{3,T}(-a,0,c),
      SVector{3,T}(0,a,c),
      SVector{3,T}(0,-a,c)) ),
    SymmetryOrbit(2,4, (a::T, c::T) -> (
      SVector{3,T}(a,a,c),
      SVector{3,T}(-a,a,c),
      SVector{3,T}(a,-a,c),
      SVector{3,T}(-a,-a,c)) ),
    SymmetryOrbit(3,8, (a::T, b::T, c::T) -> (
      SVector{3,T}(a,b,c),
      SVector{3,T}(b,a,c),
      SVector{3,T}(-a,b,c),
      SVector{3,T}(b,-a,c),
      SVector{3,T}(a,-b,c),
      SVector{3,T}(-b,a,c),
      SVector{3,T}(-a,-b,c),
      SVector{3,T}(-b,-a,c)) )]
end
