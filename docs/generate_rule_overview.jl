using FiniteElementQuadratureRules
using YAML
using Base: Filesystem
using Printf: @sprintf

_domain_symbol(::Type{D}) where {D} = Symbol(D)
_domain_name(::Type{D}) where {D} = uppercasefirst(string(Symbol(D)))
_domain_filename(::Type{D}) where {D} = lowercase(string(Symbol(D))) * ".md"

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

function _format_accuracy(x)
  value = x isa AbstractString ? parse(BigFloat, x) : BigFloat(x)
  return @sprintf("%.2e", Float64(value))
end

function _format_rule_reference(reference::AbstractString)
  key = _markdown_escape(reference)
  return "[$key](@cite)"
end

function _format_orbits(data::Dict)
  haskey(data, "orbits") || return "-"
  "[" * join(string.(data["orbits"]), ", ") * "]"
end

function _overview_header(io::IO, title::AbstractString)
  println(io, "# $(title)")
  println(io)
  println(io, "Quality markers:")
  println(io, "- `P`: positive weights")
  println(io, "- `N`: some negative weights")
  println(io, "- `I`: all quadrature points are strictly inside the reference element")
  println(io, "- `B`: some quadrature points lie on the boundary")
  println(io, "- `O`: some quadrature points lie outside the reference element")
end

function _write_domain_page(
  out_file::AbstractString,
  domain_name::AbstractString,
  domain_entries::AbstractVector;
)
  open(out_file, "w") do io
    _overview_header(io, "$(domain_name) Rules")
    println(io)
    println(io, "[Back to Database Overview](index.md)")

    if isempty(domain_entries)
      println(io)
      println(io, "_No rules found for this domain._")
      return
    end

    println(io)
    println(io, "| Degree | Points | Accuracy | Quality | Orbits | Reference |")
    println(io, "| ---: | ---: | ---: | :---: | :--- | :--- |")
    for entry in domain_entries
      println(
        io,
        "| $(entry.degree) | $(entry.points) | `",
        _markdown_escape(entry.accuracy),
        "` | $(entry.quality) | `",
        _markdown_escape(entry.orbits),
        "` | ",
        _format_rule_reference(entry.reference),
        " |",
      )
    end
  end
end

"""
    generate_rule_overview(in_dir::AbstractString, out_dir::AbstractString)

Generate a split Markdown overview of the compact quadrature-rule database, producing
one page per domain and an index page below `out_dir`.
"""
function generate_rule_overview(in_dir::AbstractString, out_dir::AbstractString)
  Filesystem.mkpath(out_dir)

  entries = Dict{Symbol, Vector{NamedTuple}}()
  for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
    entries[_domain_symbol(Domain)] = NamedTuple[]
  end

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in sort!(collect(f for f in files if endswith(f, ".yml") && !occursin(".opt.yml", f)))
      path = joinpath(root, file)
      println("Parsing '$(path)'")
      data, qr = _read_rule_file(path)
      D = Symbol(domaintype(qr))
      reference = haskey(data, "reference") ? string(data["reference"]) : "unknown"
      push!(entries[D], (
        degree = qr.degree,
        points = length(qr),
        accuracy = _format_accuracy(haskey(data, "accuracy") ? data["accuracy"] : quadratureAccuracy(qr)),
        quality = string(getQuality(qr)),
        orbits = _format_orbits(data),
        reference = reference,
        path = relpath(path, in_dir),
      ))
    end
  end

  for domain in keys(entries)
    sort!(entries[domain]; by = entry -> (entry.degree, entry.points, entry.accuracy, entry.quality, entry.orbits, entry.reference, entry.path))
  end

  generated = String[]
  for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
    D = _domain_symbol(Domain)
    out_file = joinpath(out_dir, _domain_filename(Domain))
    _write_domain_page(out_file, _domain_name(Domain), entries[D])
    push!(generated, out_file)
  end

  index_file = joinpath(out_dir, "index.md")
  open(index_file, "w") do io
    _overview_header(io, "Quadrature Rule Database")
    println(io)
    println(io, "Per-domain rule overviews:")
    println(io)
    for Domain in Base.uniontypes(FiniteElementQuadratureRules.AllDomains)
      D = _domain_symbol(Domain)
      count = length(entries[D])
      filename = _domain_filename(Domain)
      println(io, "- [", _domain_name(Domain), "](", filename, "): ", count, " rule", count == 1 ? "" : "s")
    end
  end
  push!(generated, index_file)

  generated
end
