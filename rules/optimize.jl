using FiniteElementQuadratureRules
using YAML: load_file
using Base: Filesystem
using Printf: @sprintf

const DEFAULT_ACCURACY_THRESHOLD = big"1e-20"
const DEFAULT_OPT_TOL = BigFloat(sqrt(eps(Float64)))
const DEFAULT_CLOSE_ATOL = big"1e-8"
const DEFAULT_CLOSE_RTOL = big"1e-8"

function _load_compact_rule(path::AbstractString)
  data = load_file(path)
  rule = haskey(data, "weights") ?
    CompactQuadratureRuleWithWeights(BigFloat, data) :
    CompactQuadratureRule(BigFloat, data)
  return data, rule
end

function _compact_base_rule(rule::CompactQuadratureRule)
  return rule
end

function _compact_base_rule(rule::CompactQuadratureRuleWithWeights)
  return CompactQuadratureRule(rule.domain, rule.degree, rule.orbits, rule.positions)
end

function _rebuild_rule(rule::CompactQuadratureRuleWithWeights, ocqr::CompactQuadratureRule, oqr::QuadratureRule)
  compact_weights = getCompactWeights(oqr, ocqr.orbits)
  return CompactQuadratureRuleWithWeights(rule.domain, rule.degree, ocqr.orbits, ocqr.positions, compact_weights)
end

function _rebuild_rule(::CompactQuadratureRule, ocqr::CompactQuadratureRule, ::QuadratureRule)
  return ocqr
end

function _max_componentwise_difference(a::AbstractVector, b::AbstractVector)
  length(a) == length(b) || return big(Inf)
  isempty(a) && return zero(BigFloat)
  return maximum(abs.(BigFloat.(a) .- BigFloat.(b)))
end

function _is_close(a::AbstractVector, b::AbstractVector; atol::Real=DEFAULT_CLOSE_ATOL, rtol::Real=DEFAULT_CLOSE_RTOL)
  length(a) == length(b) || return false
  isempty(a) && return true
  maxdiff = _max_componentwise_difference(a, b)
  scale = max(maximum(abs.(BigFloat.(a))), maximum(abs.(BigFloat.(b))), one(BigFloat))
  return maxdiff <= BigFloat(atol) + BigFloat(rtol) * scale
end

function _is_close(original::CompactQuadratureRule, optimized::CompactQuadratureRule; atol::Real=DEFAULT_CLOSE_ATOL, rtol::Real=DEFAULT_CLOSE_RTOL)
  original.domain == optimized.domain || return false
  original.degree == optimized.degree || return false
  original.orbits == optimized.orbits || return false
  return _is_close(original.positions, optimized.positions; atol, rtol)
end

function _is_close(original::CompactQuadratureRuleWithWeights, optimized::CompactQuadratureRuleWithWeights; atol::Real=DEFAULT_CLOSE_ATOL, rtol::Real=DEFAULT_CLOSE_RTOL)
  original.domain == optimized.domain || return false
  original.degree == optimized.degree || return false
  original.orbits == optimized.orbits || return false
  return _is_close(original.positions, optimized.positions; atol, rtol) &&
         _is_close(original.weights, optimized.weights; atol, rtol)
end

function _try_optimized_rule(rule;
                             solver=:str,
                             maxiters::Integer=120,
                             opt_tol::Real=DEFAULT_OPT_TOL,
                             show_trace::Bool=false,
                             solver_kwargs...)
  try
    base = _compact_base_rule(rule)
    ocqr = optimize(base; solver, maxiters, show_trace, solver_kwargs...)
    oqr = expand(ocqr)
    @assert testWeights(oqr; tol=opt_tol)
    @assert testQuadratureRule(oqr; tol=opt_tol)
    optimized = _rebuild_rule(rule, ocqr, oqr)
    return (
      ok = true,
      solver = solver,
      rule = optimized,
      accuracy = quadratureAccuracy(oqr),
      error = nothing,
    )
  catch err
    return (
      ok = false,
      solver = solver,
      rule = nothing,
      accuracy = big(Inf),
      error = err,
    )
  end
end

function optimize_rule_file(path::AbstractString;
                            accuracy_threshold::Real=DEFAULT_ACCURACY_THRESHOLD,
                            maxiters::Integer=120,
                            close_atol::Real=DEFAULT_CLOSE_ATOL,
                            close_rtol::Real=DEFAULT_CLOSE_RTOL,
                            opt_tol::Real=DEFAULT_OPT_TOL,
                            precision::Integer=32,
                            show_trace::Bool=false)
  data, original = _load_compact_rule(path)
  original_accuracy = haskey(data, "accuracy") ?
    parse(BigFloat, string(data["accuracy"])) :
    quadratureAccuracy(expand(original))
  reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"

  if original_accuracy <= BigFloat(accuracy_threshold)
    return (
      path = path,
      skipped = true,
      reason = :accuracy_ok,
      original_accuracy = original_accuracy,
      final_accuracy = original_accuracy,
      solver = nothing,
      action = :skip,
      output = path,
      attempts = NamedTuple[],
    )
  end

  attempted = NamedTuple[]

  str_result = _try_optimized_rule(original;
    solver = :str,
    maxiters,
    opt_tol,
    show_trace)
  attempted = vcat(attempted, [str_result])

  best = str_result
  if !str_result.ok || !(str_result.accuracy < original_accuracy)
    lm_result = _try_optimized_rule(original;
      solver = :lm,
      maxiters,
      opt_tol,
      show_trace)
    attempted = vcat(attempted, [lm_result])
    if lm_result.ok && lm_result.accuracy < min(original_accuracy, str_result.accuracy)
      best = lm_result
    end
  end

  if !best.ok || !(best.accuracy < original_accuracy)
    return (
      path = path,
      skipped = false,
      reason = :not_improved,
      original_accuracy = original_accuracy,
      final_accuracy = original_accuracy,
      solver = nothing,
      action = :skip,
      output = path,
      attempts = attempted,
    )
  end

  close = _is_close(original, best.rule; atol=close_atol, rtol=close_rtol)
  outpath = close ? path : replace(path, r"\.yml$" => ".opt.yml")
  write_file(outpath, best.rule; reference, precision)

  return (
    path = path,
    skipped = false,
    reason = :improved,
    original_accuracy = original_accuracy,
    final_accuracy = best.accuracy,
    solver = best.solver,
    action = close ? :overwrite : :optfile,
    output = outpath,
    attempts = attempted,
  )
end

function optimize_rules(dir::AbstractString=joinpath(@__DIR__, "compact");
                        accuracy_threshold::Real=DEFAULT_ACCURACY_THRESHOLD,
                        maxiters::Integer=120,
                        close_atol::Real=DEFAULT_CLOSE_ATOL,
                        close_rtol::Real=DEFAULT_CLOSE_RTOL,
                        opt_tol::Real=DEFAULT_OPT_TOL,
                        precision::Integer=32,
                        show_trace::Bool=false)
  candidates = NamedTuple[]
  for (root, _, files) in Filesystem.walkdir(dir)
    for file in sort(files)
      endswith(file, ".yml") || continue
      occursin(".opt.yml", file) && continue
      path = joinpath(root, file)
      data = load_file(path)
      accuracy = haskey(data, "accuracy") ? parse(BigFloat, string(data["accuracy"])) : big(Inf)
      accuracy > BigFloat(accuracy_threshold) || continue
      candidates = vcat(candidates, [(path = path, accuracy = accuracy)])
    end
  end
  candidates = sort(candidates; by = x -> (x.accuracy, x.path), rev = true)

  results = NamedTuple[]
  println("Candidates above accuracy threshold: $(length(candidates))")
  for candidate in candidates
    path = candidate.path
    result = optimize_rule_file(path;
      accuracy_threshold,
      maxiters,
      close_atol,
      close_rtol,
      opt_tol,
      precision,
      show_trace)
    results = vcat(results, [result])

    if result.action == :overwrite
      println(@sprintf("OVERWRITE %s  %.3e -> %.3e  [%s]",
        relpath(path, pwd()), Float64(result.original_accuracy), Float64(result.final_accuracy), string(result.solver)))
    elseif result.action == :optfile
      println(@sprintf("OPTFILE   %s  %.3e -> %.3e  [%s] -> %s",
        relpath(path, pwd()), Float64(result.original_accuracy), Float64(result.final_accuracy), string(result.solver), relpath(result.output, pwd())))
    elseif result.reason == :accuracy_ok
      println(@sprintf("SKIP      %s  %.3e", relpath(path, pwd()), Float64(result.original_accuracy)))
    else
      println(@sprintf("NOCHANGE  %s  %.3e", relpath(path, pwd()), Float64(result.original_accuracy)))
      for attempt in result.attempts
        if !attempt.ok
          println("  $(attempt.solver): ERROR $(sprint(showerror, attempt.error))")
        else
          println(@sprintf("  %s: %.3e", string(attempt.solver), Float64(attempt.accuracy)))
        end
      end
    end
  end
  return results
end

function low_accuracy_rule_files(dir::AbstractString=joinpath(@__DIR__, "compact");
                                 accuracy_threshold::Real=DEFAULT_ACCURACY_THRESHOLD)
  candidates = NamedTuple[]
  for (root, _, files) in Filesystem.walkdir(dir)
    for file in sort(files)
      endswith(file, ".yml") || continue
      occursin(".opt.yml", file) && continue
      path = joinpath(root, file)
      data = load_file(path)
      accuracy = haskey(data, "accuracy") ? parse(BigFloat, string(data["accuracy"])) : big(Inf)
      accuracy > BigFloat(accuracy_threshold) || continue
      candidates = vcat(candidates, [(path = path, accuracy = accuracy)])
    end
  end
  candidates = sort(candidates; by = x -> (x.accuracy, x.path), rev = true)
  return [candidate.path for candidate in candidates]
end

if abspath(PROGRAM_FILE) == @__FILE__
  if isempty(ARGS)
    optimize_rules()
  else
    for path in ARGS
      optimize_rule_file(path)
    end
  end
end
