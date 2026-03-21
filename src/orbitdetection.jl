"""
    detectSymmetryOrbits(domain::AbstractDomain, points::AbstractVector;
                         weights=nothing,
                         candidate_orbits=collect(eachindex(symmetryOrbits(T, domain))),
                         point_tol=sqrt(eps(T)),
                         weight_tol=point_tol,
                         canonicalize=identity)

Infer a compact orbit representation from `points`, which must already be expressed in
the coordinate system expected by the symmetry orbits of `domain` (for example,
barycentric coordinates on triangles and tetrahedra).

If `weights` are provided, only points with matching weights are grouped into the same
orbit, and one representative weight per detected orbit is returned unchanged. The
result is a named tuple with fields `orbits`, `positions`, and, when applicable,
`weights`.
"""
function detectSymmetryOrbits(domain::AbstractDomain, points::AbstractVector;
                              weights=nothing,
                              candidate_orbits=nothing,
                              point_tol::Real=sqrt(eps(Float64)),
                              weight_tol::Real=point_tol,
                              canonicalize=identity)
  isempty(points) && return isnothing(weights) ?
    (orbits=Int[], positions=Float64[]) :
    (orbits=Int[], positions=Float64[], weights=Float64[])

  T = foldl(promote_type, (eltype(p) for p in points); init=Float64)
  sos = symmetryOrbits(T, domain)
  orbit_indices = isnothing(candidate_orbits) ? collect(eachindex(sos)) : collect(candidate_orbits)

  normalized_points = [canonicalize(p) for p in points]
  normalized_weights = isnothing(weights) ? nothing : T[_transform(T, w) for w in weights]

  used = falses(length(points))
  orbit_counts = zeros(Int, length(sos))
  positions_by_orbit = [T[] for _ in eachindex(sos)]
  weights_by_orbit = [T[] for _ in eachindex(sos)]

  function point_distance(p, q)
    maximum(abs.(T.(p) .- T.(q)))
  end

  function find_matching_point!(candidate_used::BitVector, target_point, target_weight)
    for i in eachindex(normalized_points)
      candidate_used[i] && continue
      if !isnothing(normalized_weights) && abs(normalized_weights[i] - target_weight) > T(weight_tol)
        continue
      end
      if point_distance(normalized_points[i], target_point) <= T(point_tol)
        candidate_used[i] = true
        return i
      end
    end
    return nothing
  end

  while true
    idx = findfirst(!, used)
    isnothing(idx) && break

    matched = false
    for orbit_index in orbit_indices
      so = sos[orbit_index]
      params = T[compact(so, normalized_points[idx])...]
      expanded_points = [canonicalize(p) for p in collect(expand(so, params))]
      isempty(expanded_points) && continue

      candidate_used = copy(used)
      target_weight = isnothing(normalized_weights) ? zero(T) : normalized_weights[idx]
      valid = true
      for point in expanded_points
        match_idx = find_matching_point!(candidate_used, point, target_weight)
        if isnothing(match_idx)
          valid = false
          break
        end
      end

      if valid
        used = candidate_used
        orbit_counts[orbit_index] += 1
        append!(positions_by_orbit[orbit_index], params)
        if !isnothing(normalized_weights)
          push!(weights_by_orbit[orbit_index], target_weight)
        end
        matched = true
        break
      end
    end

    matched || error("Failed to match point $(normalized_points[idx]) to any candidate symmetry orbit.")
  end

  positions = reduce(vcat, positions_by_orbit; init=T[])
  if isnothing(normalized_weights)
    return (orbits=orbit_counts, positions=positions)
  end

  orbit_weights = reduce(vcat, weights_by_orbit; init=T[])
  return (orbits=orbit_counts, positions=positions, weights=orbit_weights)
end
