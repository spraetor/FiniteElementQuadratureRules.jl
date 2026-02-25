using FiniteElementQuadratureRules
using YAML: load_file
using Base: Filesystem
using Printf: @sprintf

function fix_z_coordinate()
  counter = 0
  for (root, _, files) in Filesystem.walkdir(@__DIR__)
    for file in (f for f in files if !startswith(f,"_") && endswith(f, ".yml"))
      println("Parsing '$(joinpath(root, file))'")
      data = load_file(joinpath(root, file))
      if haskey(data, "weights")
        cqr = CompactQuadratureRuleWithWeights(BigFloat, data)
      else
        cqr = CompactQuadratureRule(BigFloat, data)
      end

      sos = symmetryOrbits(ctype(cqr),cqr.domain)
      j = 1
      for (i,orbits) in enumerate(cqr.orbits)
        so = sos[i]
        n = args(so)
        for _ in 1:orbits
          point = Vector(expand(so,view(cqr.positions,j:j+n-1))[1])
          point[3] = 2*point[3]-1
          cqr.positions[j:j+n-1] .= compact(so, point)
          j = j + n
        end
      end

      open(joinpath(root, file), "w") do f
        write(f, "reference: '$(data["reference"])'\n")
        write(f, "region: $(region(cqr.domain))\n")
        write(f, "dim: $(dimension(cqr.domain))\n")
        write(f, "degree: $(cqr.degree)\n")
        write(f, "orbits: [$(cqr.orbits[1])")
        for i in 2:length(cqr.orbits)
          write(f, ", $(cqr.orbits[i])")
        end
        write(f, "]\n")
        write(f, "positions:\n")
        for p in cqr.positions
          write(f, "  - '$(@sprintf("%0.*e",50,p))'\n")
        end
      end
    end
  end
end