using Base: Filesystem
import YAML


function expandall(in_dir::AbstractString, out_dir::AbstractString)
  out_dir = Filesystem.mkpath(Filesystem.dirname(out_dir))

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      println("read $(joinpath(root, file))")
      data = YAML.load_file(joinpath(root, file))
      qr = if haskey(data, "weights")
        expand(CompactQuadratureRuleWithWeights(data))
      else
        expand(CompactQuadratureRule(data))
      end

      out_subdir = Filesystem.mkpath(root)
      out_file = joinpath(out_dir, relpath(out_subdir, in_dir), basename(file))
      println("  -> $(out_file)")
      YAML.write_file(out_file, Dict(qr, data["reference"]))
    end
  end
end
