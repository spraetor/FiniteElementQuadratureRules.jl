using YAML: load_file, write_file
using Base: Filesystem
using OteraEngine

function default_chooser(qr1::QuadratureRule{T, D, Domain}, qr2::QuadratureRule{T, D, Domain}) where {T,D,Domain}
  if qr1.degree != qr2.degree
    # compare the degree of the quadrature rule
    return qr1.degree < qr2.degree
  elseif length(qr1) != length(qr2)
    # compare the number of points
    return length(qr1) < length(qr2)
  elseif (:inside in qr1.properties) != (:inside in qr2.properties)
    # check whether a rule has only inside points
    return :inside in qr1.properties
  elseif (:positive in qr1.properties) != (:positive in qr2.properties)
    # check whether a rule has only positive weights
    return :positive in qr1.properties
  else
    return true
    # TODO: connect to BibFormatter/BibParser...
    # # compare the publication years
    # if length(qr2.bib.date.year) == 0
    #   return true
    # elseif length(qr1.bib.date.year) == 0
    #   return false
    # else
    #   year1 = parse(Int, qr1.bib.date.year)
    #   year2 = parse(Int, qr2.bib.date.year)
    #   return year1 < year2
    # end
  end
end

function generate(template::AbstractString, in_dir::AbstractString, out_dir::AbstractString;
                  filter = (qr)->true, chooser=default_chooser, kwargs...)
  out_dir = Filesystem.mkpath(out_dir)

  qrs = Dict{Symbol, Vector{Tuple{QuadratureRule, String}}}()
  for Domain in Base.uniontypes(AllDomains)
    qrs[Symbol(Domain)] = Tuple{QuadratureRule, String}[]
  end

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      println("Parsing '$(joinpath(root, file))'")
      data = load_file(joinpath(root, file))
      if haskey(data, "weights")
        cqr = CompactQuadratureRuleWithWeights(BigFloat, data)
      else
        cqr = CompactQuadratureRule(BigFloat, data)
      end
      qr = expand(cqr)
      if filter(qr)
        reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
        push!(qrs[Symbol(domaintype(qr))], (qr, reference))
      end
    end
  end

  println("Sorting the quadrature rules")
  for domain in map(Symbol, Base.uniontypes(AllDomains))
    sort!(qrs[domain]; lt=(a,b)->chooser(a[1], b[1]))
    selected = Tuple{QuadratureRule, String}[]
    seen_degrees = Set{Int}()
    for rule in qrs[domain]
      degree = rule[1].degree
      if degree in seen_degrees
        continue
      end
      push!(selected, rule)
      push!(seen_degrees, degree)
    end
    qrs[domain] = selected
  end

  maxdegree = Dict{Symbol, Int}()
  for domain in map(Symbol, Base.uniontypes(AllDomains))
    maxdegree[domain] = maximum(rule[1].degree for rule in qrs[domain]; init=0)
  end

  tmpl = Template(
    template;
    config=Dict(
      "autoescape" => false,
      "autospace" => false,
      "lstrip_blocks" => false,
      "trim_blocks" => false
    )
  )

  println("Generating the output files")
  for domain in Base.uniontypes(AllDomains)
    D = Symbol(domain)
    rules = Vector{Dict{String, Any}}()
    for (qr, reference) in qrs[D]
      rule = Dict{String, Any}(Dict(qr; reference=reference, kwargs...))
      push!(rules, rule)
    end
    data = Dict(
      :domain => uppercasefirst(string(D)),
      :date => Libc.strftime("%Y-%m-%d", time()),
      :dim => dimension(domain),
      :region => uppercasefirst(region(domain)),
      :maxdegree => maxdegree[D],
      :rules => rules
    )

    open(joinpath(out_dir, string(D) * ".hh"), "w") do f
      write(f, tmpl(init=data))
    end
  end

end



function expandall(in_dir::AbstractString, out_dir::AbstractString)
  out_dir = Filesystem.mkpath(Filesystem.dirname(out_dir))

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      println("read $(joinpath(root, file))")
      data = YAML.load_file(joinpath(root, file))
      out_root = joinpath(out_dir, relpath(root, in_dir))
      out_file = joinpath(out_root, file)
      if haskey(data, "weights")
        cqr = CompactQuadratureRuleWithWeights(Float64, data)
      else
        cqr = CompactQuadratureRule(Float64, data)
      end
      qr = expand(cqr)
      if !isnothing(qr)
        mkpath(out_root)
        YAML.write_file(out_file, Dict(qr, data["reference"]))
      end
    end
  end
end
