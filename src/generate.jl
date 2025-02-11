using YAML: load_file
using Base: Filesystem

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
    # compare the publication years
    if length(qr2.bib.date.year) == 0
      return true
    elseif length(qr1.bib.date.year) == 0
      return false
    else
      year1 = parse(Int, qr1.bib.date.year)
      year2 = parse(Int, qr2.bib.date.year)
      return year1 < year2
    end
  end
end

function generate(template::AbstractString, in_dir::AbstractString, out_dir::AbstractString;
                  filter = (qr)->true, chooser=default_chooser, kwargs...)
  out_dir = Filesystem.mkpath(Filesystem.dirname(out_dir))

  qrs = Dict{Symbol, Vector{QuadratureRule}}()
  for Domain in Base.uniontypes(AllDomains)
    qrs[Symbol(Domain)] = QuadratureRule{String,dimension(Domain),Domain}[]
  end
  println(typeof(qrs))

  for (root, _, files) in Filesystem.walkdir(in_dir)
    for file in (f for f in files if endswith(f, ".yml"))
      data = load_file(joinpath(root, file))
      dim = data["dim"]
      qr = QuadratureRule{String}(dim, data["region"], data; kwargs...)
      if filter(qr)
        push!(qrs[Symbol(domaintype(qr))], qr)
      end
    end
  end

  for domain in map(Symbol, Base.uniontypes(AllDomains))
    sort!(qrs[domain]; lt=chooser)
  end

  println(qrs)

  maxdegree = Dict{Symbol, Int}()
  for domain in map(Symbol, Base.uniontypes(AllDomains))
    maxdegree[domain] = maximum(qr.degree for qr in qrs[domain]; init=0)
  end
  println(maxdegree)
end