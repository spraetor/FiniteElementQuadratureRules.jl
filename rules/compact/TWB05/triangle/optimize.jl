using FiniteElementQuadratureRules
using YAML: load_file
using Base: Filesystem
using Printf: @sprintf

const TWB05_REFERENCE = "TWB05"
const OPT_TOL = sqrt(eps(Float64))

function optimize_twb05_rule(file_in::AbstractString, file_out::AbstractString=file_in;
                             maxiters::Integer=120, solver=:lm, solver_kwargs...)
  data = load_file(file_in)
  cqr = CompactQuadratureRuleWithWeights(BigFloat, data)
  base = CompactQuadratureRule(cqr.domain, cqr.degree, cqr.orbits, cqr.positions)
  ocqr = optimize(base; show_trace = true, maxiters = maxiters,
                  solver = solver, solver_kwargs...)
  oqr = expand(ocqr)
  @assert testWeights(oqr; tol = OPT_TOL)
  @assert testQuadratureRule(oqr; tol = OPT_TOL)

  compact_weights = getCompactWeights(oqr, ocqr.orbits)
  optimized = CompactQuadratureRuleWithWeights(cqr.domain, cqr.degree, ocqr.orbits, ocqr.positions, compact_weights)
  write_file(file_out, optimized; reference = data["reference"], precision = 32)
  return file_out
end

function optimize_twb05_rules(dir_in::AbstractString=@__DIR__, dir_out::AbstractString=dir_in;
                              maxiters::Integer=120, solver=:lm, solver_kwargs...)
  mkpath(dir_out)
  outputs = String[]
  failures = Pair{String,Any}[]
  for file in sort(filter(f -> endswith(f, ".yml"), readdir(dir_in; join=true)))
    outfile = joinpath(dir_out, basename(file))
    tmpfile = outfile * ".tmp"
    try
      println("Optimizing '$(basename(file))'")
      push!(outputs, optimize_twb05_rule(file, tmpfile;
                                         maxiters, solver, solver_kwargs...))
      Filesystem.mv(tmpfile, outfile; force=true)
    catch err
      isfile(tmpfile) && rm(tmpfile; force=true)
      push!(failures, file => err)
      println("Failed to optimize '$(basename(file))': $(sprint(showerror, err))")
    end
  end
  if !isempty(failures)
    println("\nOptimization failures:")
    for (file, err) in failures
      println(@sprintf("  %s: %s", basename(file), sprint(showerror, err)))
    end
  end
  return outputs
end
