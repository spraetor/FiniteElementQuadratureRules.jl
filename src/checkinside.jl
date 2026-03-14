"""
    checkInside(ref::ReferenceElement, x::AbstractVector, tol::Real)
    checkInside(ref::ReferenceElement, x::AbstractVector{T})

Check whether a point `x` lies inside the reference element `ref` or on its boundary.
This check is performed with a given tolerance `tol`. If no tolerance is given, use
the tolerance `tol=eps(T)`.
"""
function checkInside end

"""
    checkStrictlyInside(ref::ReferenceElement, x::AbstractVector, tol::Real)
    checkStrictlyInside(ref::ReferenceElement, x::AbstractVector{T})

Check whether a point `x` lies strictly inside the reference element `ref`.
This check is performed with a given tolerance `tol`. If no tolerance is given, use
the tolerance `tol=eps(T)`.
"""
function checkStrictlyInside end


function checkInside(ref::ReferenceElement, x::AbstractVector{T}) where {T<:Real}
  tol::T = (T <: AbstractFloat) ? eps(T) : zero(T)
  checkInside(ref,x,tol)
end

checkInside(::ReferenceElement{Point}, ::AbstractVector, ::Real) = true

function checkInside(ref::ReferenceElement{Line}, x::AbstractVector, tol::Real)
  coordinates(ref, 1)[1]-tol <= x[1] <= coordinates(ref, 2)[1]+tol
end

function checkInside(ref::ReferenceElement{<:AbstractSimplex}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(ref, x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ)
end

function checkInside(ref::ReferenceElement{<:AbstractCube}, x::AbstractVector, tol::Real)
  all(coordinates(ref, 1) .- tol .<= x .<= coordinates(ref)[end] .+ tol)
end

function checkInside(ref::ReferenceElement{Prism}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(ref, x)
  all(-tol <= λᵢ <= 1+tol for λᵢ in λ[1:3]) && (-1-tol <= λ[4] <= 1+tol)
end

# TODO: generalize implementation to arbitrary reference pyramids
function checkInside(::ReferenceElement{Pyramid}, x::AbstractVector, tol::Real)
  z = x[3]
  s = (1 - z)/2
  -1-tol <= z <= 1+tol && abs(x[1]) <= s + tol && abs(x[2]) <= s + tol
end



function checkStrictlyInside(ref::ReferenceElement, x::AbstractVector{T}) where {T<:Real}
  tol::T = (T <: AbstractFloat) ? eps(T) : zero(T)
  checkStrictlyInside(ref,x,tol)
end

checkStrictlyInside(::ReferenceElement{Point}, ::AbstractVector, ::Real) = false

function checkStrictlyInside(ref::ReferenceElement{Line}, x::AbstractVector, tol::Real)
  coordinates(ref, 1)[1]+tol < x[1] < coordinates(ref, 2)[1]-tol
end

function checkStrictlyInside(ref::ReferenceElement{<:AbstractSimplex}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(ref, x)
  all(tol < λᵢ < 1-tol for λᵢ in λ)
end

function checkStrictlyInside(ref::ReferenceElement{<:AbstractCube}, x::AbstractVector, tol::Real)
  all(coordinates(ref, 1) .+ tol .< x .< coordinates(ref)[end] .- tol)
end

function checkStrictlyInside(ref::ReferenceElement{Prism}, x::AbstractVector, tol::Real)
  λ = barycentricCoordinates(ref, x)
  all(tol < λᵢ < 1-tol for λᵢ in λ[1:3]) && (-1+tol < λ[4] < 1-tol)
end

# TODO: generalize implementation to arbitrary reference pyramids
function checkStrictlyInside(::ReferenceElement{Pyramid}, x::AbstractVector, tol::Real)
  z = x[3]
  s = (1 - z)/2
  -1+tol < z < 1-tol && abs(x[1]) < s - tol && abs(x[2]) < s - tol
end
