using StaticArrays: SVector

function duneReferenceElement(::Point)
  ReferenceElement(Point, [SVector{0,Int}()], Vector{Int}[], 1)
end

function duneReferenceElement(::Line)
  ReferenceElement(Line, SVector{1,Int}[[0], [1]], [[1], [2]], 1)
end

function duneReferenceElement(::Triangle)
  ReferenceElement(Triangle, SVector{2,Int}[[0,0], [1,0], [0,1]],
    [[1,2], [1,3], [2,3]], 1//2)
end

function duneReferenceElement(::Quadrilateral)
  ReferenceElement(Quadrilateral, SVector{2,Int}[[0,0], [1,0], [0,1], [1,1]],
    [[1,3], [2,4], [1,2], [3,4]], 1)
end

function duneReferenceElement(::Tetrahedron)
  ReferenceElement(Tetrahedron, SVector{3,Int}[[0,0,0], [1,0,0], [0,1,0], [0,0,1]],
    [[1,2,3], [1,2,4], [1,3,4], [2,3,4]], 1//6)
end

function duneReferenceElement(::Hexahedron)
  ReferenceElement(Hexahedron, SVector{3,Int}[[0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1], [1,0,1], [0,1,1], [1,1,1]],
    [[1,3,5,7], [2,4,6,8], [1,2,5,6], [3,4,7,8], [1,2,3,4], [6,7,8,9]], 1)
end

function duneReferenceElement(::Prism)
  ReferenceElement(Prism, SVector{3,Int}[[0,0,0], [1,0,0], [0,1,0], [0,0,1], [1,0,1], [0,1,1]],
    [[1,2,4,5], [1,3,4,6], [2,3,5,6], [1,2,3], [4,5,6]], 1//2)
end

function duneReferenceElement(::Pyramid)
  ReferenceElement(Pyramid,SVector{3,Int}[[0,0,0], [1,0,0], [0,1,0], [1,1,0], [0,0,1]],
    [[1,2,3,4], [1,3,5], [2,4,5], [1,2,5], [3,4,5]], 1//3)
end
