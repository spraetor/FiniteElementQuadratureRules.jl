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
    sort!(entries[domain]; by = entry -> (entry.degree, entry.points, entry.quality, entry.orbits, entry.reference, entry.path))
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
