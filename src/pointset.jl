struct PointSet{Point}
  points::Vector{Point}
end

using StaticArrays: MVector, @SVector
function EquidistantPointSet(::Type{ct}, Ω::Triangle, n) where {ct<:Real}
  h = one(ct)/n
  points = [ zeros(ct,3) for _ in 1:(n+2)*(n+1)/2 ]
  k = 1
  for i in 0:n
    for j = 0:n-i
      points[k] .= (1-(i+j)*h, j*h, i*h)
      k = k+1
    end
  end

  points
end

function LagrangeBasis(Ω::Triangle, n)
  points = EquidistantPointSet(Rational, Ω, n)
  indices = Vector{Int}[ Int.(p.*n) for p in points ]
  basis0 = Function[ λ -> prod(prod(λ[j] - p//n for p in 0:I[j]-1; init=1) for j in eachindex(I); init=1) for I in indices ]
  Function[ λ -> b(λ)//b(λᵢ) for (b,λᵢ) in zip(basis0, points) ]
end

using Symbolics: @variables
function SymbolicLagrangeBasis(Ω::Triangle, n)
  λ = @variables λ₁ λ₂ λ₃
  [ b(λ) for b in LagrangeBasis(Ω,n) ]
end

using Symbolics: coeff, degree
function integrate(Ω::Triangle, f)
  λ = @variables λ₁ λ₂ λ₃
  c = [ coeff(f,λᵢ) for λᵢ in (λ₁,λ₂,λ₃) ]
end