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
function makeCompactQuadratureRule(::Type{T}, dim::Integer, region::AbstractString, data::Dict) where T<:Real
  dom = domain(dim,region)
  degree = Int(data["degree"])
  orbits = Int[ o for o in data["orbits"] ]
  if haskey(data, "positions")
    positions = isnothing(data["positions"]) ? T[] : T[ _parse(T,p) for p in data["positions"] ]
  elseif haskey(data, "arguments")
    positions = isnothing(data["arguments"]) ? T[] : T[ _parse(T,p) for p in data["arguments"] ]
  else
    return nothing
  end

  CompactQuadratureRule{typeof(dom),T}(dom,degree,orbits,positions)
end

ctype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω


# expand a compact rule into a quadrature rule
function expand(cqr::CompactQuadratureRule{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  if length(cqr.orbits) > length(sos)
    error("Number of orbits incompatible with available symmetric orbits.")
    return nothing
  end

  j = 1
  Point = typeof(sos[1].expand()[1])
  points = Point[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = so.args
    for _ in 1:orbits
      push!(points, so.expand(cqr.positions[j:j+n-1]...)...)
      j = j + n
    end
  end

  # maybe compute the weights here directly
  coords = transformcoordinates(cqr.domain,points)
  weights = getWeights(T,cqr.domain,cqr.degree,coords,cqr.orbits)
  QuadratureRule(cqr.domain, cqr.degree, coords, weights)
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
function makeCompactQuadratureRuleWithWeights(::Type{T}, dim::Integer, region::AbstractString, data::Dict) where T
  dom = domain(dim,region)
  degree = Int(data["degree"])
  orbits = Int[ o for o in data["orbits"] ]
  positions = isnothing(data["positions"]) ? T[] : T[ _parse(T,p) for p in data["positions"] ]
  weights = isnothing(data["weights"]) ? T[] : T[ _parse(T,w) for w in data["weights"] ]

  CompactQuadratureRuleWithWeights(dom,degree,orbits,positions,weights)
end

ctype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = T
domaintype(::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real} = Ω

# expand a compact rule with weights into a full quadrature rule
function expand(cqr::CompactQuadratureRuleWithWeights{Ω,T}) where {Ω<:AbstractDomain,T<:Real}
  sos = symmetryOrbits(T,cqr.domain)
  if length(cqr.orbits) > length(sos)
    return nothing
  end
  j = 1
  k = 1
  Point = typeof(sos[1].expand()[1])
  points = Point[]
  weights = T[]
  for (i,orbits) in enumerate(cqr.orbits)
    so = sos[i]
    n = so.args
    for _ in 1:orbits
      push!(points, so.expand(cqr.positions[j:j+n-1]...)...)
      append!(weights, fill(cqr.weights[k],so.size))
      j = j + n
      k = k + 1
    end
  end

  @assert length(points) == length(weights)
  QuadratureRule(cqr.domain, cqr.degree,
    transformcoordinates(cqr.domain,points), transformweights(cqr.domain, weights))
end



function expandall(in_dir::AbstractString, out_dir::AbstractString)
  out_dir = Filesystem.mkpath(Filesystem.dirname(out_dir))

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      println("read $(joinpath(root, file))")
      data = load_file(joinpath(root, file))
      out_root = joinpath(out_dir, relpath(root, in_dir))
      out_file = joinpath(out_root, file)
      dim = data["dim"]
      if haskey(data, "weights")
        cqr = makeCompactQuadratureRuleWithWeights(Float64, dim, data["region"], data)
      else
        cqr = makeCompactQuadratureRule(Float64, dim, data["region"], data)
      end
      qr = expand(cqr)
      if !isnothing(qr)
        mkpath(out_root)
        write_file(out_file, qr, data)
      end
    end
  end
end