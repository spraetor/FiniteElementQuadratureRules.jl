using StaticArray: SVector

struct ReferenceElement{D, Ω<:AbstractDomain}
  coordinates::Vector{SVector{D,Int}}
  facets::Vector{Vector{Int}}
end

# Define standard reference elements

function ReferenceElement(::Point)
  ReferenceElement{0,Point}([[]], [])
end

function ReferenceElement(::Line)
  ReferenceElement{1,Line}([[0], [1]], [[1], [2]])
end

function ReferenceElement(::Triangle)
  ReferenceElement{2,Triangle}([[0,0], [1,0], [0,1]],
    [[1,2], [1,3], [2,3]])
end

function ReferenceElement(::Quadrilateral)
  ReferenceElement{2,Quadrilateral}([[0,0], [1,0], [0,1], [1,1]],
    [[1,3], [2,4], [1,2], [3,4]])
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

dimension(::ReferenceElement{D,Ω}) where {D, Ω<:AbstractDomain} = D
domaintype(::ReferenceElement{D,Ω}) where {D, Ω<:AbstractDomain} = Ω

volume(::Type{T}, ::ReferenceElement{2,Triangle}) = T(1)
volume(::Type{T}, ::ReferenceElement{1,Line}) = T(1)
volume(::Type{T}, ::ReferenceElement{2,Triangle}) = T(1//2)
volume(::Type{T}, ::ReferenceElement{2,Quadrilateral}) = T(1)
volume(::Type{T}, ::ReferenceElement{3,Tetrahedron}) = T(1//6)
volume(::Type{T}, ::ReferenceElement{3,Hexahedron}) = T(1)
volume(::Type{T}, ::ReferenceElement{3,Prism}) = T(1//2)
volume(::Type{T}, ::ReferenceElement{3,Pyramid}) = T(1//3)

volume(ref::ReferenceElement{D,Ω}) where {D,Ω<:AbstractDomain} = volume(Rational,ref)