using StaticArrays: SVector

"""
    ReferenceElement{Ω,Point}

A reference representation of a given domain Ω in terms of corner coordinates
and a list of facets connecting the corners.
"""
struct ReferenceElement{Ω<:AbstractDomain, Point<:AbstractVector}
  coordinates::Vector{Point}
  facets::Vector{Vector{Int}}
end

dimension(::ReferenceElement{Ω,P}) where {Ω,P} = dimension(Ω)
domaintype(::ReferenceElement{Ω,P}) where {Ω,P} = Ω
domain(::ReferenceElement{Ω,P}) where {Ω,P} = Ω()

function ReferenceElement(::Point)
  ReferenceElement{Point,SVector{0,Int}}([[]], [])
end

function ReferenceElement(::Line)
  ReferenceElement{Line,SVector{1,Int}}([[-1], [1]], [[1], [2]])
end

function ReferenceElement(::Triangle)
  ReferenceElement{Triangle,SVector{2,Int}}([[-1,-1], [1,-1], [-1,1]],
    [[1,2], [1,3], [2,3]])
end

function ReferenceElement(::Quadrilateral)
  ReferenceElement{Quadrilateral,SVector{2,Int}}([[-1,-1], [1,-1], [-1,1], [1,1]],
    [[1,3], [2,4], [1,2], [3,4]])
end

function ReferenceElement(::Tetrahedron)
  ReferenceElement{Tetrahedron,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1]],
    [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])
end

function ReferenceElement(::Hexahedron)
  ReferenceElement{Hexahedron,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [-1,-1,1], [1,-1,1], [-1,1,1], [1,1,1]],
    [[1,3,5,7], [2,4,6,8], [1,2,5,6], [3,4,7,8], [1,2,3,4], [5,6,7,8]])
end

function ReferenceElement(::Prism)
  ReferenceElement{Prism,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1], [1,-1,1], [-1,1,1]],
    [[1,2,4,5], [1,3,4,6], [2,3,5,6], [1,2,3], [4,5,6]])
end

function ReferenceElement(::Pyramid)
  ReferenceElement{Pyramid,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [0,0,1]],
    [[1,2,3,4], [1,3,5], [2,4,5], [1,2,5], [3,4,5]])
end

volume(::Type{T}, ::ReferenceElement{Point,P}) where {T<:Real,P} = T(1)
volume(::Type{T}, ::ReferenceElement{Line,P}) where {T<:Real,P} = T(2)
volume(::Type{T}, ::ReferenceElement{Triangle,P}) where {T<:Real,P} = T(2)
volume(::Type{T}, ::ReferenceElement{Quadrilateral,P}) where {T<:Real,P} = T(4)
volume(::Type{T}, ::ReferenceElement{Tetrahedron,P}) where {T<:Real,P} = T(4//3)
volume(::Type{T}, ::ReferenceElement{Hexahedron,P}) where {T<:Real,P} = T(8)
volume(::Type{T}, ::ReferenceElement{Prism,P}) where {T<:Real,P} = T(4)
volume(::Type{T}, ::ReferenceElement{Pyramid,P}) where {T<:Real,P} = T(8//3)

volume(ref::ReferenceElement{Ω,P}) where {Ω,P} = volume(Rational,ref)


import Base: position

"""
    position(ref::ReferenceElement{Ω,P}, i::Integer, c::Integer)

Compute the center of the `i`th sub-entity of codimension `c` of the reference
element `ref`.
"""
function position(ref::ReferenceElement{Ω,P}, i::Integer, c::Integer) where {Ω<:AbstractDomain, P}
  if dimension(Ω) == 2
    _position2d(ref,i,c)
  elseif dimension(Ω) == 3
    _position3d(ref,i,c)
  else
    error("Not implemented!")
  end
end

function position(ref::ReferenceElement{Point,P}, i::Integer, c::Integer) where P
  @assert c == 0 && i == 1
  return ref.coordinates[1]
end

function position(ref::ReferenceElement{Line,P}, i::Integer, c::Integer) where P
  @assert c <= 1
  if c == 0
    @assert i == 1
    return (ref.coordinates[1] + ref.coordinates[2])/2
  elseif c == 1
    @assert i <= 2
    return ref.coordinates[i]
  end
end


function _position2d(ref::ReferenceElement, i::Integer, c::Integer)
  @assert dimension(ref) == 2
  @assert c <= 2
  if c == 0
    @assert i == 1
    return sum(c for c in ref.coordinates)/length(ref.coordinates)
  elseif c == 1
    @assert i <= length(ref.facets)
    facet = ref.facets[i]
    return sum(ref.coordinates[i] for i in facet)/length(facet)
  elseif c == 2
    @assert i <= length(ref.coordinates)
    return ref.coordinates[i]
  end
end

function _position3d(ref::ReferenceElement, i::Integer, c::Integer)
  @assert dimension(ref) == 3
  @assert c <= 3
  if c == 0
    @assert i == 1
    return sum(c for c in ref.coordinates)/length(ref.coordinates)
  elseif c == 1
    @assert i <= length(ref.facets)
    facet = ref.facets[i]
    return sum(ref.coordinates[i] for i in facet)/length(facet)
  elseif c == 2
    error("Not implemented.")
  elseif c == 3
    @assert i <= length(ref.coordinates)
    return ref.coordinates[i]
  end
end


"""
    checkInside(ref::ReferenceElement, x::AbstractVector, tol::Real)

Check whether a point `x` lies inside the reference element `ref` or on its boundary.
This check is performed with a given tolerance `tol`.
"""
function checkInside end

"""
    checkInside(ref::ReferenceElement, x::AbstractVector{T})

Check whether a point `x` lies inside the reference element `ref` or on its boundary.
This check is performed with a tolerance `tol=eps(T)`.
"""
function checkInside(ref::ReferenceElement, x::AbstractVector{T}) where {T<:Real}
  tol::T = (T <: AbstractFloat) ? eps(T) : zero(T)
  checkInside(ref,x,tol)
end

checkInside(::ReferenceElement{Point,P}, ::AbstractVector, ::Real) where P = true

function checkInside(ref::ReferenceElement{Line,P}, x::AbstractVector, tol::Real) where P
  ref.coordinates[1][1]-tol <= x[1] <= ref.coordinates[2][1]+tol
end

function checkInside(::ReferenceElement{Ω,P}, x::AbstractVector, tol::Real) where {Ω<:AbstractSimplex,P}
  λ = barycentricCoordinates(Ω(), x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ)
end

function checkInside(ref::ReferenceElement{Ω,P}, x::AbstractVector, tol::Real) where {Ω<:AbstractCube,P}
  all(ref.coordinates[1] .- tol .<= x .<= ref.coordinates[end] .+ tol)
end

function checkInside(::ReferenceElement{Prism,P}, x::AbstractVector, tol::Real) where P
  λ = barycentricCoordinates(Prism(), x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ[1:3]) && (-1-tol <= λ[4] <= 1+tol)
end

function checkInside(::ReferenceElement{Pyramid,P}, x::AbstractVector, tol::Real) where P
  z = x[3]
  s = (1 - z)/2
  -1-tol <= z <= 1+tol && abs(x[1]) <= s + tol && abs(x[2]) <= s + tol
end


"""
    checkStrictlyInside(ref::ReferenceElement, x::AbstractVector, tol::Real)

Check whether a point `x` lies strictly inside the reference element `ref`.
This check is performed with a given tolerance `tol`.
"""
function checkStrictlyInside end

"""
    checkStrictlyInside(ref::ReferenceElement, x::AbstractVector{T})

Check whether a point `x` lies strictly inside the reference element `ref`.
This check is performed with a tolerance `tol=eps(T)`.
"""
function checkStrictlyInside(ref::ReferenceElement, x::AbstractVector{T}) where {T<:Real}
  tol::T = (T <: AbstractFloat) ? eps(T) : zero(T)
  checkStrictlyInside(ref,x,tol)
end


checkStrictlyInside(::ReferenceElement{Point,P}, ::AbstractVector, ::Real) where P = false

function checkStrictlyInside(ref::ReferenceElement{Line,P}, x::AbstractVector, tol::Real) where P
  ref.coordinates[1][1]+tol < x[1] < ref.coordinates[2][1]-tol
end

function checkStrictlyInside(::ReferenceElement{Ω,P}, x::AbstractVector, tol::Real) where {Ω<:AbstractSimplex,P}
  λ = barycentricCoordinates(Ω(), x)
  all(tol < λᵢ < 1-tol for λᵢ in λ)
end

function checkStrictlyInside(ref::ReferenceElement{Ω,P}, x::AbstractVector, tol::Real) where {Ω<:AbstractCube,P}
  all(ref.coordinates[1] .+ tol .< x .< ref.coordinates[end] .- tol)
end

function checkStrictlyInside(::ReferenceElement{Prism,P}, x::AbstractVector, tol::Real) where P
  λ = barycentricCoordinates(Prism(), x)
  all(tol < λᵢ < 1-tol for λᵢ in λ[1:3]) && (-1+tol < λ[4] < 1-tol)
end

function checkStrictlyInside(::ReferenceElement{Pyramid,P}, x::AbstractVector, tol::Real) where P
  z = x[3]
  s = (1 - z)/2
  -1+tol < z < 1-tol && abs(x[1]) < s - tol && abs(x[2]) < s - tol
end
