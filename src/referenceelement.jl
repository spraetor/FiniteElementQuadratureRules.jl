using StaticArrays: SVector

"""
  ReferenceElement{D,Ω}
"""
struct ReferenceElement{D, Ω<:AbstractDomain}
  coordinates::Vector{SVector{D,Int}}
  facets::Vector{Vector{Int}}
end

dimension(::ReferenceElement{D,Ω}) where {D, Ω<:AbstractDomain} = D
domaintype(::ReferenceElement{D,Ω}) where {D, Ω<:AbstractDomain} = Ω

# Define standard reference elements

function ReferenceElement(::Point)
  ReferenceElement{0,Point}([[]], [])
end

import Base: position
function position(ref::ReferenceElement{0,Point}, i::Integer, c::Integer)
  @assert c == 0 && i == 1
  return ref.coordinates[1]
end

function ReferenceElement(::Line)
  ReferenceElement{1,Line}([[0], [1]], [[1], [2]])
end

function position(ref::ReferenceElement{1,Line}, i::Integer, c::Integer)
  @assert c <= 1
  if c == 0
    @assert i == 1
    return (ref.coordinates[1] + ref.coordinates[2])/2
  elseif c == 1
    @assert i <= 2
    return ref.coordinates[i]
  end
end

function ReferenceElement(::Triangle)
  ReferenceElement{2,Triangle}([[0,0], [1,0], [0,1]],
    [[1,2], [1,3], [2,3]])
end

function ReferenceElement(::Quadrilateral)
  ReferenceElement{2,Quadrilateral}([[0,0], [1,0], [0,1], [1,1]],
    [[1,3], [2,4], [1,2], [3,4]])
end

function position(ref::ReferenceElement{2,<:AbstractDomain}, i::Integer, c::Integer)
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

function ReferenceElement(::Tetrahedron)
  ReferenceElement{3,Tetrahedron}([[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    [[1,2,3], [1,2,4], [1,3,4], [2,3,4]])
end

function ReferenceElement(::Hexahedron)
  ReferenceElement{3,Hexahedron}([[0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1], [1,0,1], [0,1,1], [1,1,1]],
    [[1,3,5,7], [2,4,6,8], [1,2,5,6], [3,4,7,8], [1,2,3,4], [6,7,8,9]])
end

function ReferenceElement(::Prism)
  ReferenceElement{3,Prism}([[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,0,1], [0,1,1]],
    [[1,2,4,5], [1,3,4,6], [2,3,5,6], [1,2,3], [4,5,6]])
end

function ReferenceElement(::Pyramid)
  ReferenceElement{3,Pyramid}([[0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1]],
    [[1,2,3,4], [1,3,5], [2,4,5], [1,2,5], [3,4,5]])
end

function position(ref::ReferenceElement{3,<:AbstractDomain}, i::Integer, c::Integer)
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

volume(::Type{T}, ::ReferenceElement{0,Point}) where T<:Real = T(1)
volume(::Type{T}, ::ReferenceElement{1,Line}) where T<:Real = T(1)
volume(::Type{T}, ::ReferenceElement{2,Triangle}) where T<:Real = T(1//2)
volume(::Type{T}, ::ReferenceElement{2,Quadrilateral}) where T<:Real = T(1)
volume(::Type{T}, ::ReferenceElement{3,Tetrahedron}) where T<:Real = T(1//6)
volume(::Type{T}, ::ReferenceElement{3,Hexahedron}) where T<:Real = T(1)
volume(::Type{T}, ::ReferenceElement{3,Prism}) where T<:Real = T(1//2)
volume(::Type{T}, ::ReferenceElement{3,Pyramid}) where T<:Real = T(1//3)

volume(ref::ReferenceElement{D,Ω}) where {D,Ω<:AbstractDomain} = volume(Rational,ref)