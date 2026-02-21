using StaticArrays: SVector

"""
  ReferenceElement{D,Ω}
"""
struct ReferenceElement{D, Ω<:AbstractDomain, Point<:AbstractVector}
  coordinates::Vector{Point}
  facets::Vector{Vector{Int}}
end

dimension(::ReferenceElement{D,Ω,P}) where {D,Ω,P} = D
domaintype(::ReferenceElement{D,Ω,P}) where {D,Ω,P} = Ω

function ReferenceElement(::Point)
  ReferenceElement{0,Point,SVector{0,Int}}([[]], [])
end

function ReferenceElement(::Line)
  ReferenceElement{1,Line,SVector{1,Int}}([[-1], [1]], [[1], [2]])
end

function ReferenceElement(::Triangle)
  ReferenceElement{2,Triangle,SVector{2,Int}}([[-1,-1], [1,-1], [-1,1]],
    [[1,2], [1,3], [2,3]])
end

function ReferenceElement(::Quadrilateral)
  ReferenceElement{2,Quadrilateral,SVector{2,Int}}([[-1,-1], [1,-1], [-1,1], [1,1]],
    [[1,3], [2,4], [1,2], [3,4]])
end

function ReferenceElement(::Tetrahedron)
  ReferenceElement{3,Tetrahedron,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1]],
    [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])
end

function ReferenceElement(::Hexahedron)
  ReferenceElement{3,Hexahedron,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [-1,-1,1], [1,-1,1], [-1,1,1], [1,1,1]],
    [[1,3,5,7], [2,4,6,8], [1,2,5,6], [3,4,7,8], [1,2,3,4], [5,6,7,8]])
end

function ReferenceElement(::Prism)
  ReferenceElement{3,Prism,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1], [1,-1,1], [-1,1,1]],
    [[1,2,4,5], [1,3,4,6], [2,3,5,6], [1,2,3], [4,5,6]])
end

function ReferenceElement(::Pyramid)
  ReferenceElement{3,Pyramid,SVector{3,Int}}([[-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [0,0,1]],
    [[1,2,3,4], [1,3,5], [2,4,5], [1,2,5], [3,4,5]])
end

volume(::Type{T}, ::ReferenceElement{0,Point,P}) where {T<:Real,P} = T(1)
volume(::Type{T}, ::ReferenceElement{1,Line,P}) where {T<:Real,P} = T(2)
volume(::Type{T}, ::ReferenceElement{2,Triangle,P}) where {T<:Real,P} = T(2)
volume(::Type{T}, ::ReferenceElement{2,Quadrilateral,P}) where {T<:Real,P} = T(4)
volume(::Type{T}, ::ReferenceElement{3,Tetrahedron,P}) where {T<:Real,P} = T(4//3)
volume(::Type{T}, ::ReferenceElement{3,Hexahedron,P}) where {T<:Real,P} = T(8)
volume(::Type{T}, ::ReferenceElement{3,Prism,P}) where {T<:Real,P} = T(4)
volume(::Type{T}, ::ReferenceElement{3,Pyramid,P}) where {T<:Real,P} = T(8//3)

volume(ref::ReferenceElement{D,Ω,P}) where {D,Ω,P} = volume(Rational,ref)


import Base: position
function position(ref::ReferenceElement{0,Point,P}, i::Integer, c::Integer) where P
  @assert c == 0 && i == 1
  return ref.coordinates[1]
end

function position(ref::ReferenceElement{1,Line,P}, i::Integer, c::Integer) where P
  @assert c <= 1
  if c == 0
    @assert i == 1
    return (ref.coordinates[1] + ref.coordinates[2])/2
  elseif c == 1
    @assert i <= 2
    return ref.coordinates[i]
  end
end

function position(ref::ReferenceElement{2,<:AbstractDomain,P}, i::Integer, c::Integer) where P
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

function position(ref::ReferenceElement{3,<:AbstractDomain,P}, i::Integer, c::Integer) where P
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


checkInside(::ReferenceElement{0,Point,P}, ::AbstractVector, ::Real) where P = true

function checkInside(ref::ReferenceElement{1,Line,P}, x::AbstractVector, tol::Real) where P
  ref.coordinates[1][1]-tol <= x[1] <= ref.coordinates[2][1]+tol
end

function checkInside(ref::ReferenceElement{dim,Ω,P}, x::AbstractVector, tol::Real) where {dim,Ω<:AbstractSimplex,P}
  λ = barycentricCoordinates(Ω(), x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ)
end

function checkInside(ref::ReferenceElement{dim,Ω,P}, x::AbstractVector, tol::Real) where {dim,Ω<:AbstractCube,P}
  all(ref.coordinates[1] .- tol .<= x .<= ref.coordinates[end] .+ tol)
end

function checkInside(::ReferenceElement{3,Prism,P}, x::AbstractVector, tol::Real) where P
  λ = barycentricCoordinates(Prism(), x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ[1:3]) && (-1-tol <= λ[4] <= 1+tol)
end

function checkInside(::ReferenceElement{3,Pyramid,P}, x::AbstractVector, tol::Real) where P
  z = x[3]
  s = (1 - z)/2
  -1-tol <= z <= 1+tol && abs(x[1]) <= s + tol && abs(x[2]) <= s + tol
end

# Check whether a point x is inside the domain of the reference element `ref` or on the boundary
function checkInside(ref::ReferenceElement, x::AbstractVector{T}) where {T<:Real}
  tol::T = (T <: AbstractFloat) ? eps(T) : zero(T)
  checkInside(ref,x,tol)
end


checkStrictlyInside(::ReferenceElement{0,Point,P}, ::AbstractVector, ::Real) where P = false

function checkStrictlyInside(ref::ReferenceElement{1,Line,P}, x::AbstractVector, tol::Real) where P
  ref.coordinates[1][1]+tol < x[1] < ref.coordinates[2][1]-tol
end

function checkStrictlyInside(ref::ReferenceElement{dim,Ω,P}, x::AbstractVector, tol::Real) where {dim,Ω<:AbstractSimplex,P}
  λ = barycentricCoordinates(Ω(), x)
  all(tol < λᵢ < 1-tol for λᵢ in λ)
end

function checkStrictlyInside(ref::ReferenceElement{dim,Ω,P}, x::AbstractVector, tol::Real) where {dim,Ω<:AbstractCube,P}
  all(ref.coordinates[1] .+ tol .< x .< ref.coordinates[end] .- tol)
end

function checkStrictlyInside(::ReferenceElement{3,Prism,P}, x::AbstractVector, tol::Real) where P
  λ = barycentricCoordinates(Prism(), x)
  all(tol < λᵢ < 1-tol for λᵢ in λ[1:3]) && (-1+tol < λ[4] < 1-tol)
end

function checkStrictlyInside(::ReferenceElement{3,Pyramid,P}, x::AbstractVector, tol::Real) where P
  z = x[3]
  s = (1 - z)/2
  -1+tol < z < 1-tol && abs(x[1]) < s - tol && abs(x[2]) < s - tol
end

# Check whether a point x is strictly inside the domain of the reference element `ref`
function checkStrictlyInside(ref::ReferenceElement, x::AbstractVector{T}) where {T<:Real}
  tol::T = (T <: AbstractFloat) ? eps(T) : zero(T)
  checkStrictlyInside(ref,x,tol)
end