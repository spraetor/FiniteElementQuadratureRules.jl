"""
  CompactQuadratureRule{T,Ω}
"""
struct CompactQuadratureRule{Ω<:AbstractDomain, T<:Real}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  positions::Vector{T}
end

# Create a quadrature rule from a Dict of strings and string arrays
function (CompactQuadratureRule{T})(dim::Integer, region::AbstractString, data::Dict) where T
  D::Val = Val(dim)

  domain = domain(D,region)
  degree = Int(data["degree"])
  orbits = Int[ o for o in data["orbits"] ]
  positions = T[ _parse(T,p) for p in data["positions"] ]

  CompactQuadratureRule(domain,degree,orbits,positions)
end

ctype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω

# expand a compact rule into a quadrature rule
function expand(cqr::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  j = 1
  Point = typeof(sos[1].orbit()[1])
  points = Point[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = so.args
    for _ in 1:orbits
      push!(points, so.orbit(cqr.positions[j:j+n-1]...)...)
      j = j + n
    end
  end

  # maybe compute the weights here directly
  weights = getWeights(T,cqr.domain,cqr.degree,points,cqr.orbits)
  QuadratureRule(cqr.domain, cqr.degree, points, weights)
end

"""
  CompactQuadratureRuleWithWeights{T,Ω}

Compact rule that stores also weights.
"""
struct CompactQuadratureRuleWithWeights{Ω<:AbstractDomain, T<:Real}
  domain::Ω
  degree::Int
  orbits::Vector{Int}
  positions::Vector{T}
  weights::Vector{T}
end

# Create a quadrature rule from a Dict of strings and string arrays
function (CompactQuadratureRuleWithWeights{T})(dim::Integer, region::AbstractString, data::Dict) where T
  D::Val = Val(dim)

  domain = domain(D,region)
  degree = Int(data["degree"])
  orbits = Int[ o for o in data["orbits"] ]
  positions = T[ _parse(T,p) for p in data["positions"] ]
  weights = T[ _parse(T,w) for w in data["weights"] ]

  CompactQuadratureRuleWithWeights(domain,degree,orbits,positions,weights)
end

ctype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω

# expand a compact rule with weights into a full quadrature rule
function expand(cqr::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  j = 1
  k = 1
  Point = typeof(sos[1].orbit()[1])
  points = Point[]
  weights = T[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = so.args
    for _ in 1:orbits
      push!(points, so.orbit(cqr.positions[j:j+n-1]...)...)
      append!(weights, fill(cqr.weights[k],so.size))
      j = j + n
      k = k + 1
    end
  end

  @assert length(points) == length(weights)
  QuadratureRule(cqr.domain, cqr.degree, points, weights)
end



function expandall(in_dir::AbstractString, out_dir::AbstractString)
  out_dir = Filesystem.mkpath(Filesystem.dirname(out_dir))

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      data = load_file(joinpath(root, file))
      out_file = joinpath(out_dir, relpath(root, in_dir), file)
      dim = data["dim"]
      if haskey(data, "weights")
        cqr = CompactQuadratureRuleWithWeights{String}(dim, data["region"], data)
        write_file(out_file, expand(cqr), data)
      else
        cqr = CompactQuadratureRule{String}(dim, data["region"], data)
        write_file(out_file, expand(cqr), data)
      end
    end
  end
end