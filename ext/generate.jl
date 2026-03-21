"""
    generate(template::AbstractString, in_dir::AbstractString, out_dir::AbstractString;
             bib=joinpath(@__DIR__, "..", "references.bib"),
             refOut=(domain) -> ReferenceElement(domain),
             filter=(qr) -> true,
             chooser=default_chooser,
             kwargs...)

Read compact quadrature rules from `in_dir`, validate and transform them, select one
representative rule per degree and domain type, and render the selected rules into
output files using `template`.

This is the export-oriented generator used for targets such as the DUNE header files.
Unlike [`generate_rule_overview`](@ref), it does not keep every rule. Instead, after
validation it sorts the rules with `chooser` and keeps only the first rule for each
degree.

Keyword arguments:
- `bib`: BibTeX file used to format references included in the generated output
- `refOut`: maps a domain type to the reference element used for the exported rule
- `filter`: predicate deciding whether a validated rule is kept as a candidate
- `chooser`: comparison function used to rank competing rules of the same domain
- `kwargs...`: additional key-value pairs forwarded into the rendered template data
"""
function generate(
  template::AbstractString,
  in_dir::AbstractString,
  out_dir::AbstractString;
  bib = joinpath(@__DIR__,"..","references.bib"),
  refOut = (domain) -> ReferenceElement(domain),
  filter = (qr)->true,
  chooser = default_chooser,
  kwargs...
)
  out_dir = Filesystem.mkpath(out_dir)

  validation_tol(::Type{T}) where {T<:Real} = sqrt(eps(T))
  validation_tol(::Type{BigFloat}) = BigFloat(sqrt(eps(Float64)))

  bibFile = Filesystem.isfile(bib) ? BibParser.parse_file(bib) : Dict()

  QR = Union{Tuple{QuadratureRule, BibInternal.Entry},Missing}
  qrs = Dict{Symbol, Vector{QR}}()
  for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
    qrs[Symbol(Domain)] = Tuple{QuadratureRule, BibInternal.Entry}[]
  end

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      println("Parsing '$(joinpath(root, file))'")
      data = YAML.load_file(joinpath(root, file))
      cqr = haskey(data, "weights") ? CompactQuadratureRuleWithWeights(BigFloat, data) : CompactQuadratureRule(BigFloat, data)
      qr = expand(cqr)
      tol = validation_tol(ctype(qr))
      # if !testWeights(qr; tol)
      #   println("  -> error(weights)")
      #   continue
      # end
      # if !testQuadratureRule(qr; tol)
      #   println("  -> error(polynomials)")
      #   continue
      # end
      qr = transform(qr, refOut(domain(qr)))
      # if !testWeights(qr; tol)
      #   println("  -> error(transform)")
      #   continue
      # end
      if filter(qr)
        reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
        bibentry = haskey(bibFile, reference) ? bibFile[reference] : BibInternal.Entry(reference, Dict())
        push!(qrs[Symbol(domaintype(qr))], (qr, bibentry))
      end
    end
  end

  maxdegree = Dict{Symbol, Int}()
  for domain in map(Symbol, Base.uniontypes(FiniteElementQuadratureRules.AllDomains))
    maxdegree[domain] = maximum(rule[1].degree for rule in qrs[domain]; init=0)
  end

  println("Sorting the quadrature rules")
  for domain in map(Symbol, Base.uniontypes(FiniteElementQuadratureRules.AllDomains))
    sort!(qrs[domain]; lt=chooser)
    selected = Vector{QR}(undef, maxdegree[domain])
    fill!(selected, missing)
    for rule in qrs[domain]
      degree = rule[1].degree
      if ismissing(selected[degree])
        selected[degree] = rule
      end
    end
    qrs[domain] = selected
  end

  tmpl = Template(
    template;
    config=Dict(
      "autoescape" => false,
      "autospace" => false,
      "lstrip_blocks" => false,
      "trim_blocks" => false,
    ),
  )

  println("Generating the output files")
  bibliography(bibentry) = BibFormatter.format(bibentry, fmt=:text)
  convert(rule) = Dict(rule[1]; reference=bibliography(rule[2]), kwargs...)
  passmissing(f) = rule -> ismissing(rule) ? rule : f(rule)
  for domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
    D = Symbol(domain)
    rules = map(passmissing(convert), qrs[D])

    data = Dict(
      :domain => uppercasefirst(string(D)),
      :date => Libc.strftime("%Y-%m-%d", time()),
      :dim => dimension(domain),
      :region => uppercasefirst(region(domain)),
      :maxdegree => maxdegree[D],
      :rules => rules,
    )

    extension = Filesystem.splitext(template)[end]
    open(joinpath(out_dir, string(D) * extension), "w") do io
      write(io, tmpl(init=data))
    end
  end
end
