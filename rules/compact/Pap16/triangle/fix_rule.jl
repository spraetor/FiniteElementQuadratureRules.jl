using FiniteElementQuadratureRules
using Base: Filesystem
using Double64
import YAML

function fix_rule(file, fileOut=file)
  data = YAML.load_file(file)
  if haskey(data, "weights")
    cqr = CompactQuadratureRuleWithWeights(Double64, data)
  else
    cqr = CompactQuadratureRule(Double64, data)
  end

  ocqr = optimize(cqr)
  write_file(fileOut, ocqr)
end