using YAML: load_file, write_file
@testset "README Examples" begin
  rules_root = joinpath(@__DIR__, "..", "rules")

  @testset "Expand Compact Rule" begin
    data = load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(BigFloat, data)
    qr = expand(cqr)

    @test qr.degree == 4
    @test length(qr) == 6

    mktempdir() do tmp
      out = joinpath(tmp, "4-6.yml")
      write_file(out, Dict(qr; reference=data["reference"]))
      @test isfile(out)
      exported = load_file(out)
      @test Int(exported["degree"]) == 4
      @test length(exported["coordinates"]) == 6
    end
  end

  @testset "Transform Rule To Dune Convention" begin
    data = load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    qr = expand(CompactQuadratureRule(Float64, data))

    ref_dune = duneReferenceElement(domain(qr))
    qr_dune = transform(qr, ref_dune)
    @test qr_dune.degree == qr.degree
    @test length(qr_dune) == length(qr)
  end

  @testset "Generate Dune Files From Template" begin
    mktempdir() do tmp
      out_dir = joinpath(tmp, "dune")
      generate(
        joinpath(@__DIR__, "..", "dune.templ.hh"),
        joinpath(rules_root, "compact", "Gat88"),
        out_dir;
        precision=80
      )

      triangle_hh = joinpath(out_dir, "triangle.hh")
      tetrahedron_hh = joinpath(out_dir, "tetrahedron.hh")
      @test isfile(triangle_hh)
      @test isfile(tetrahedron_hh)

      triangle_text = read(triangle_hh, String)
      @test occursin("QuadratureRule", triangle_text)
      @test occursin("highest_order", triangle_text)
    end
  end

  @testset "Optimize Rule" begin
    data = load_file(joinpath(rules_root, "compact", "CCGV22", "triangle", "4-6.yml"))
    cqr = CompactQuadratureRule(BigFloat, data)
    oqr = optimize(cqr)
    @test oqr.degree == cqr.degree
    @test length(oqr.positions) == length(cqr.positions)
  end
end
