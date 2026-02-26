using FiniteElementQuadratureRules
using Base: Filesystem
using Quadmath: Float128
using DoubleFloats: Double64
import YAML

function fix_rule(::Type{T}, file, fileOut=file) where T
  data = YAML.load_file(file)
  if haskey(data, "weights")
    cqr = CompactQuadratureRuleWithWeights(T, data)
  else
    cqr = CompactQuadratureRule(T, data)
  end
  # cqr.positions .+= (2*rand(T,length(cqr.positions)) .- T(1))./100 .* cqr.positions
  ocqr = optimize(cqr, show_trace=true)
  qr = expand(ocqr)
  testWeights(qr)
  testQuadratureRule(qr)
  write_file(fileOut, ocqr, reference=data["reference"])
end

function fix_rule(file, fileOut=file)
  fix_rule(Double64, file, fileOut)
end

function fix_rules()
  for (root, _, files) in Filesystem.walkdir(@__DIR__)
    for file in (f for f in files if startswith(f,"rot") && endswith(f, ".yml"))
      if Filesystem.isfile(joinpath(root, file*"_"))
        continue
      end
      println("Optimizing '$(file)'")
      fix_rule(Double64,joinpath(root, file),joinpath(root, file*"_"))
    end
  end
end

function copy_rules()
  for (root, _, files) in Filesystem.walkdir(@__DIR__)
    for file in (f for f in files if startswith(f,"rot") && endswith(f, ".yml_"))
      println("Moving '$(file)'")
      Filesystem.mv(joinpath(root,file), joinpath(root,replace(file,"yml_"=>"yml")))
    end
  end
end
