struct PointSet{Point}
  points::Vector{Point}
end

function EquidistantPointSet{ct}(Ω<:Triangle, n) where {ct<:Real}
  h = one(ct)/n
  Point = StaticVector{3,ct}
  points = Point[ Point(0,0,0) for _ in 1:(n+2)*(n+1)/2 ]
  k = 1
  for i in 0:n
    for j = 0:i
      points[k] .= (i*h, j*h, 1-(i+j)*h)
      k = k+1
    end
  end

  PointSet{Point}(points)
end