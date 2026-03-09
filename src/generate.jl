using Base: Filesystem
using OteraEngine
import BibParser
import BibInternal
import BibFormatter
import YAML

function default_chooser(qr1_entry, qr2_entry)
  qr1, qr2 = qr1_entry[1], qr2_entry[1]
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
    entry1, entry2 = qr1_entry[2], qr2_entry[2]
    # compare the publication years
    if length(entry2.date.year) == 0
      return true
    elseif length(entry1.date.year) == 0
      return false
    else
      year1 = parse(Int, entry1.date.year)
      year2 = parse(Int, entry2.date.year)
      return year1 < year2
    end
  end
end

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
  # Decimal YAML data is parsed as BigFloat, but its source precision is typically
  # around Float64/printed-decimal accuracy rather than full BigFloat machine epsilon.
  validation_tol(::Type{BigFloat}) = BigFloat(sqrt(eps(Float64)))

  # Read a bibliography file for including a full reference
  bibFile = Filesystem.isfile(bib) ? BibParser.parse_file(bib) : Dict()

  QR = Union{Tuple{QuadratureRule, BibInternal.Entry},Missing}
  qrs = Dict{Symbol, Vector{QR}}()
  for Domain in Base.uniontypes(AllDomains)
    qrs[Symbol(Domain)] = Tuple{QuadratureRule, String}[]
  end

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      println("Parsing '$(joinpath(root, file))'")
      data = YAML.load_file(joinpath(root, file))
      if haskey(data, "weights")
        cqr = CompactQuadratureRuleWithWeights(BigFloat, data)
      else
        cqr = CompactQuadratureRule(BigFloat, data)
      end
      qr = expand(cqr)
      tol = validation_tol(ctype(qr))
      if !testWeights(qr; tol)
        println("  -> error(weights)")
        continue
      end
      if !testQuadratureRule(qr; tol)
        println("  -> error(polynomials)")
        continue
      end
      qr = transform(qr, refOut(domain(qr)))
      if !testWeights(qr; tol)
        println("  -> error(transform)")
        continue
      end
      if filter(qr)
        reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
        bibentry = haskey(bibFile, reference) ? bibFile[reference] : BibInternal.Entry(reference, Dict())
        push!(qrs[Symbol(domaintype(qr))], (qr, bibentry))
      end
    end
  end

  maxdegree = Dict{Symbol, Int}()
  for domain in map(Symbol, Base.uniontypes(AllDomains))
    maxdegree[domain] = maximum(rule[1].degree for rule in qrs[domain]; init=0)
  end

  println("Sorting the quadrature rules")
  for domain in map(Symbol, Base.uniontypes(AllDomains))
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
      "trim_blocks" => false
    )
  )

  println("Generating the output files")
  bibliography(bibentry) = BibFormatter.format(bibentry, fmt=:text)
  convert(r) = Dict(r[1]; reference=bibliography(r[2]), kwargs...)
  passmissing(f) = r->ismissing(r) ? r : f(r)
  for domain in Base.uniontypes(AllDomains)
    D = Symbol(domain)
    rules = map(passmissing(convert), qrs[D])

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


function _read_rule_file(path::AbstractString)
  data = YAML.load_file(path)
  qr = if haskey(data, "coordinates")
    QuadratureRule(BigFloat, data)
  elseif haskey(data, "weights")
    expand(CompactQuadratureRuleWithWeights(BigFloat, data))
  else
    expand(CompactQuadratureRule(BigFloat, data))
  end
  return data, qr
end


_markdown_escape(s::AbstractString) = replace(s, "|" => "\\|", "\n" => " ")
_markdown_escape(x) = _markdown_escape(string(x))


function _format_rule_reference(reference::AbstractString, bibfile::AbstractDict)
  label = "`$(_markdown_escape(reference))`"
  haskey(bibfile, reference) || return label
  return "[$label](#$(_markdown_escape(reference)))"
end


function _format_orbits(data::Dict)
  haskey(data, "orbits") || return "-"
  "[" * join(string.(data["orbits"]), ", ") * "]"
end


"""
    generate_rule_overview(in_dir::AbstractString, out_file::AbstractString;
                           bib=joinpath(@__DIR__, "..", "references.bib"),
                           style=:abbrv)

Read all YAML rules found below `in_dir` and generate a Markdown overview table at
`out_file`. In contrast to [`generate`](@ref), this function keeps all rules and does
not select only one representative per degree.

Rules are grouped by domain type and listed with:
- degree
- number of quadrature points after expansion
- quality marker
- orbit signature (for compact rules)
- reference key

The output also includes a bibliography section. Reference keys in the tables link to
their formatted bibliography entry in the same Markdown file.
"""
function generate_rule_overview(
  in_dir::AbstractString,
  out_file::AbstractString;
  bib = joinpath(@__DIR__,"..","references.bib"),
  style::Symbol = :abbrv,
)
  bibfile = Filesystem.isfile(bib) ? BibParser.parse_file(bib) : Dict()

  entries = Dict{Symbol, Vector{NamedTuple}}()
  for Domain in Base.uniontypes(AllDomains)
    entries[Symbol(Domain)] = NamedTuple[]
  end

  references = Set{String}()
  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in sort!(collect(f for f in files if endswith(f, ".yml")))
      path = joinpath(root, file)
      println("Parsing '$(path)'")
      data, qr = _read_rule_file(path)
      reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
      push!(references, reference)
      push!(entries[Symbol(domaintype(qr))], (
        degree = qr.degree,
        points = length(qr),
        quality = string(getQuality(qr)),
        orbits = _format_orbits(data),
        reference = reference,
        path = relpath(path, in_dir),
      ))
    end
  end

  for domain in keys(entries)
    sort!(entries[domain]; by = e -> (e.degree, e.points, e.quality, e.orbits, e.reference, e.path))
  end

  Filesystem.mkpath(Filesystem.dirname(out_file))
  open(out_file, "w") do io
    println(io, "# Quadrature Rule Overview")
    println(io)
    println(io, "This page lists all quadrature rules found in `$(in_dir)`.")
    println(io)
    println(io, "Quality markers:")
    println(io, "- `P`: positive weights")
    println(io, "- `N`: some negative weights")
    println(io, "- `I`: all quadrature points are strictly inside the reference element")
    println(io, "- `B`: some quadrature points lie on the boundary")
    println(io, "- `O`: some quadrature points lie outside the reference element")

    for Domain in Base.uniontypes(AllDomains)
      D = Symbol(Domain)
      isempty(entries[D]) && continue

      println(io)
      println(io, "## $(uppercasefirst(string(D)))")
      println(io)
      println(io, "| Degree | Points | Quality | Orbits | Reference |")
      println(io, "| ---: | ---: | :---: | :--- | :--- |")
      for entry in entries[D]
        println(
          io,
          "| $(entry.degree) | $(entry.points) | $(entry.quality) | `",
          _markdown_escape(entry.orbits),
          "` | ",
          _format_rule_reference(entry.reference, bibfile),
          " |",
        )
      end
    end

    if !isempty(references)
      println(io)
      println(io, "## References")
      for reference in sort!(collect(references))
        println(io)
        println(io, "### $(reference)")
        println(io)
        if haskey(bibfile, reference)
          println(io, BibFormatter.format(bibfile[reference]; style=style, fmt=:md))
        else
          println(io, "_Missing bibliography entry._")
        end
      end
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
