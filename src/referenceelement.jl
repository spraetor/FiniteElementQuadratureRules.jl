using StaticArrays: SVector

"""
    ReferenceElement{Ω,T,P}

A reference representation of a given domain Ω in terms of corner coordinates
and a list of facets connecting the corners.
"""
struct ReferenceElement{Ω<:AbstractDomain, T<:Real, P<:AbstractVector}
  coordinates::Vector{P}
  facets::Vector{Vector{Int}}
  volume::T

  function ReferenceElement{Ω,T,P}(coordinates::AbstractVector{Q}, facets, volume) where {Ω<:AbstractDomain, T<:Real, P<:AbstractVector, Q<:AbstractVector}
    if P != Q
      new([P(x) for x in coordinates],facets,T(volume))
    else
      new(coordinates,facets,T(volume))
    end
  end

  # converting constructor
  function ReferenceElement{Ω,T,P}(ref::ReferenceElement{Ω,S,Q}) where {Ω<:AbstractDomain, T<:Real, P<:AbstractVector, S<:Real, Q<:AbstractVector}
    ReferenceElement{Ω,T,P}(ref.coordinates, ref.facets, ref.volume)
  end
end

ReferenceElement(::Type{Ω}, coordinates::AbstractVector{P}, facets::Vector, volume::T) where {Ω<:AbstractDomain, P<:AbstractVector, T<:Real} = ReferenceElement{Ω,T,P}(coordinates,facets,volume)

function ReferenceElement(::Point)
  ReferenceElement(Point, [SVector{0,Int}()], Vector{Int}[], 1)
end

function ReferenceElement(::Line)
  ReferenceElement(Line, SVector{1,Int}[[-1], [1]], [[1], [2]], 2)
end

function ReferenceElement(::Triangle)
  ReferenceElement(Triangle, SVector{2,Int}[[-1,-1], [1,-1], [-1,1]],
    [[1,2], [1,3], [2,3]], 2)
end

function ReferenceElement(::Quadrilateral)
  ReferenceElement(Quadrilateral, SVector{2,Int}[[-1,-1], [1,-1], [-1,1], [1,1]],
    [[1,3], [2,4], [1,2], [3,4]], 4)
end

function ReferenceElement(::Tetrahedron)
  ReferenceElement(Tetrahedron, SVector{3,Int}[[-1,-1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1]],
    [[1,2,3], [1,2,4], [1,3,4], [2,3,4]], 4//3)
end

function ReferenceElement(::Hexahedron)
  ReferenceElement(Hexahedron, SVector{3,Int}[[-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [-1,-1,1], [1,-1,1], [-1,1,1], [1,1,1]],
    [[1,3,5,7], [2,4,6,8], [1,2,5,6], [3,4,7,8], [1,2,3,4], [5,6,7,8]], 8)
end

function ReferenceElement(::Prism)
  ReferenceElement(Prism, SVector{3,Int}[[-1,-1,-1], [1,-1,-1], [-1,1,-1], [-1,-1,1], [1,-1,1], [-1,1,1]],
    [[1,2,4,5], [1,3,4,6], [2,3,5,6], [1,2,3], [4,5,6]], 4)
end

function ReferenceElement(::Pyramid)
  ReferenceElement(Pyramid, SVector{3,Int}[[-1,-1,-1], [1,-1,-1], [-1,1,-1], [1,1,-1], [0,0,1]],
    [[1,2,3,4], [1,3,5], [2,4,5], [1,2,5], [3,4,5]], 8//3)
end


dimension(::ReferenceElement{Ω}) where {Ω} = dimension(Ω)
domaintype(::ReferenceElement{Ω}) where {Ω} = Ω
domain(::ReferenceElement{Ω}) where {Ω} = Ω()
ctype(::ReferenceElement{Ω,T,P}) where {Ω,T,P} = eltype(P)

coordinates(ref::ReferenceElement) = ref.coordinates
coordinates(ref::ReferenceElement, i::Integer) = ref.coordinates[i]
facets(ref::ReferenceElement) = ref.facets
facets(ref::ReferenceElement, i::Integer) = ref.facets[i]
volume(ref::ReferenceElement) = ref.volume

import Base: size

_numEdges(::Tetrahedron) = 6
_numEdges(::Pyramid) = 8
_numEdges(::Prism) = 9
_numEdges(::Hexahedron) = 12

function size(ref::ReferenceElement, codim::Integer)
  @assert codim <= dimension(ref)
  if codim == 0
    return 1
  elseif codim == 1
    return facets(domain(ref))
  elseif codim == dimension(ref)
    return vertices(domain(ref))
  else
    return _numEdges(domain(ref))
  end
end

# number of subentities of codimension cc of subentity (i,c)
function size(ref::ReferenceElement, i::Integer, c::Integer, cc::Integer)
  @assert c <= dimension(ref)
  @assert cc <= dimension(ref)
  if cc < c
    return 0
  end

  if c == 0 && i == 1
    return size(ref, cc)
  elseif c == 1
    if cc == dimension(ref)
      return length(facets(ref,i))
    elseif cc == dimension(ref)-1
      return dimension(ref)-cc == 1 ? 1 : length(facets(ref,i))
    else
      return 0
    end
  else
    return 1+cc-c
  end
end

import Base: position

"""
    position(ref::ReferenceElement, i::Integer, codim::Integer)

Compute the center of the `i`th sub-entity of codimension `codim` of the reference
element `ref`.
"""
function position(ref::ReferenceElement, i::Integer, codim::Integer)
  if dimension(ref) == 2
    _position2d(ref,i,codim)
  elseif dimension(ref) == 3
    _position3d(ref,i,codim)
  else
    error("Not implemented!")
  end
end

function position(ref::ReferenceElement{Point}, i::Integer, codim::Integer)
  @assert codim == 0 && i == 1
  return coordinates(ref, 1)
end

function position(ref::ReferenceElement{Line}, i::Integer, codim::Integer)
  @assert codim <= 1
  if codim == 0
    @assert i == 1
    return (coordinates(ref, 1) + coordinates(ref, 2))/2
  elseif codim == 1
    @assert i <= 2
    return coordinates(ref, i)
  end
end


function _position2d(ref::ReferenceElement, i::Integer, codim::Integer)
  @assert dimension(ref) == 2
  @assert codim <= 2
  if codim == 0
    @assert i == 1
    return sum(c for c in coordinates(ref))/length(coordinates(ref))
  elseif codim == 1
    @assert i <= length(facets(ref))
    facet = facets(ref, i)
    return sum(coordinates(ref, i) for i in facet)/length(facet)
  elseif codim == 2
    @assert i <= length(coordinates(ref))
    return coordinates(ref, i)
  end
end

function _position3d(ref::ReferenceElement, i::Integer, codim::Integer)
  @assert dimension(ref) == 3
  @assert codim <= 3
  if codim == 0
    @assert i == 1
    return sum(c for c in coordinates(ref))/length(coordinates(ref))
  elseif codim == 1
    @assert i <= length(facets(ref))
    facet = facets(ref, i)
    return sum(coordinates(ref, i) for i in facet)/length(facet)
  elseif codim == 2
    error("Not implemented.")
  elseif codim == 3
    @assert i <= length(coordinates(ref))
    return coordinates(ref, i)
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

checkInside(::ReferenceElement{Point}, ::AbstractVector, ::Real) = true

function checkInside(ref::ReferenceElement{Line}, x::AbstractVector, tol::Real)
  coordinates(ref, 1)[1]-tol <= x[1] <= coordinates(ref, 2)[1]+tol
end

# TODO: generalize implementation to arbitrary reference simplices
function checkInside(ref::ReferenceElement{<:AbstractSimplex}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(domain(ref), x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ)
end

function checkInside(ref::ReferenceElement{<:AbstractCube}, x::AbstractVector, tol::Real)
  all(coordinates(ref, 1) .- tol .<= x .<= coordinates(ref)[end] .+ tol)
end

# TODO: generalize implementation to arbitrary reference prisms
function checkInside(::ReferenceElement{Prism}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(Prism(), x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ[1:3]) && (-1-tol <= λ[4] <= 1+tol)
end

# TODO: generalize implementation to arbitrary reference pyramids
function checkInside(::ReferenceElement{Pyramid}, x::AbstractVector, tol::Real)
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


checkStrictlyInside(::ReferenceElement{Point}, ::AbstractVector, ::Real) = false

function checkStrictlyInside(ref::ReferenceElement{Line}, x::AbstractVector, tol::Real)
  coordinates(ref, 1)[1]+tol < x[1] < coordinates(ref, 2)[1]-tol
end

# TODO: generalize implementation to arbitrary reference simplices
function checkStrictlyInside(ref::ReferenceElement{<:AbstractSimplex}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(domain(ref), x)
  all(tol < λᵢ < 1-tol for λᵢ in λ)
end

function checkStrictlyInside(ref::ReferenceElement{<:AbstractCube}, x::AbstractVector, tol::Real)
  all(coordinates(ref, 1) .+ tol .< x .< coordinates(ref)[end] .- tol)
end

# TODO: generalize implementation to arbitrary reference prisms
function checkStrictlyInside(::ReferenceElement{Prism}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(Prism(), x)
  all(tol < λᵢ < 1-tol for λᵢ in λ[1:3]) && (-1+tol < λ[4] < 1-tol)
end

# TODO: generalize implementation to arbitrary reference pyramids
function checkStrictlyInside(::ReferenceElement{Pyramid}, x::AbstractVector, tol::Real)
  z = x[3]
  s = (1 - z)/2
  -1+tol < z < 1-tol && abs(x[1]) < s - tol && abs(x[2]) < s - tol
end
