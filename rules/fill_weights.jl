using FiniteElementQuadratureRules
using YAML: load_file
using Base: Filesystem
using Printf: @sprintf

const DEFAULT_WEIGHT_TOL = BigFloat(sqrt(eps(Float64)))

function missing_weight_rule_files(dir::AbstractString=joinpath(@__DIR__, "compact"))
  files = String[]
  for (root, _, names) in Filesystem.walkdir(dir)
    for name in sort(names)
      endswith(name, ".yml") || continue
      occursin(".opt.yml", name) && continue
      path = joinpath(root, name)
      data = load_file(path)
      haskey(data, "weights") && continue
      push!(files, path)
    end
  end
  return files
end

function fill_rule_weights(path::AbstractString;
                           precision::Integer=32,
                           tol::Real=DEFAULT_WEIGHT_TOL)
  data = load_file(path)
  haskey(data, "weights") && return (
    path = path,
    action = :skip,
    reason = :already_weighted,
  )

  cqr = CompactQuadratureRule(BigFloat, data)
  qr = expand(cqr)
  @assert testWeights(qr; tol=tol)
  @assert testQuadratureRule(qr; tol=tol)

  compact_weights = getCompactWeights(qr, cqr.orbits)
  weighted = CompactQuadratureRuleWithWeights(cqr.domain, cqr.degree, cqr.orbits, cqr.positions, compact_weights)

  reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
  write_file(path, weighted; reference, precision)

  return (
    path = path,
    action = :overwrite,
    reason = :filled,
  )
end

function fill_missing_weights(dir::AbstractString=joinpath(@__DIR__, "compact");
                              precision::Integer=32,
                              tol::Real=DEFAULT_WEIGHT_TOL)
  results = NamedTuple[]
  files = missing_weight_rule_files(dir)
  println("Rules missing weights: $(length(files))")
  for path in files
    result = fill_rule_weights(path; precision, tol)
    results = vcat(results, [result])
    if result.action == :overwrite
      println("FILLED    $(relpath(path, pwd()))")
    else
      println("SKIP      $(relpath(path, pwd()))")
    end
  end
  return results
end

if abspath(PROGRAM_FILE) == @__FILE__
  if isempty(ARGS)
    fill_missing_weights()
  else
    for path in ARGS
      fill_rule_weights(path)
    end
  end
end
