using FiniteElementQuadratureRules
using YAML: load_file
using Base: Filesystem

const RESERVED_YAML_FIELDS = Set([
  "reference", "region", "dim", "degree", "quality", "accuracy",
  "properties", "coordinates", "weights", "orbits", "positions",
])

function _extra_yaml_fields(data::AbstractDict)
  Dict{String,Any}(string(k) => v for (k, v) in pairs(data) if !(string(k) in RESERVED_YAML_FIELDS))
end

function compact_rule_files(dir::AbstractString=joinpath(@__DIR__, "compact"))
  files = String[]
  for (root, _, names) in Filesystem.walkdir(dir)
    for name in sort(names)
      endswith(name, ".yml") || continue
      path = joinpath(root, name)
      data = load_file(path)
      haskey(data, "coordinates") && continue
      push!(files, path)
    end
  end
  files
end

function update_rule_accuracy(path::AbstractString; precision::Integer=32)
  data = load_file(path)
  haskey(data, "coordinates") && return (
    path = path,
    action = :skip,
    reason = :expanded_rule,
  )

  rule = haskey(data, "weights") ?
    CompactQuadratureRuleWithWeights(BigFloat, data) :
    CompactQuadratureRule(BigFloat, data)

  reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
  extra_fields = _extra_yaml_fields(data)
  write_file(path, rule; reference, precision, extra_fields)

  return (
    path = path,
    action = :overwrite,
    reason = :updated,
  )
end

function update_accuracy(dir::AbstractString=joinpath(@__DIR__, "compact"); precision::Integer=32)
  files = compact_rule_files(dir)
  results = NamedTuple[]
  println("Compact rules: $(length(files))")
  for path in files
    result = update_rule_accuracy(path; precision)
    push!(results, result)
    if result.action == :overwrite
      println("UPDATED   $(relpath(path, pwd()))")
    else
      println("SKIP      $(relpath(path, pwd()))")
    end
  end
  results
end

if abspath(PROGRAM_FILE) == @__FILE__
  if isempty(ARGS)
    update_accuracy()
  else
    for path in ARGS
      update_rule_accuracy(path)
    end
  end
end
