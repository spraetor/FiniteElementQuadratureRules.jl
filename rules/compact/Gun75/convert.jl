using FiniteElementQuadratureRules
using StaticArrays: SVector

function gun75(file)
  # Q = sum_{i=1}^4 wt_i * f(t_i,t_i) + sum_{j=1}^5 wxy_j * (f(x_j,y_j) + f(y_j,x_j))

  t = [0.607800961271,0.275713993942,0.572041162029e-1,0.802234574054]
  wt = [-0.141499310464e-2,0.594598907982e-1,0.214817467226e-1,0.300946300607e-4]

  xy = [0.290343325778 0.627085788906e-1;
        0.623760718495 0.398844552695e-1;
        0.865863389281 0.655834126334e-1;
        0.540734125182 0.233967966580;
        0.637988448494 0.343693043712]
  wxy = [0.469862540914e-1,0.346857035218e-1,0.284354430188e-1,0.700271893727e-1,0.300870404721e-1]

  refIn = ReferenceElement(Triangle, SVector{2,Rational}[[0,0], [1,0], [0,1]],
    [[1,2], [1,3], [2,3]], 1//2)

  orbits = [0,0,0,0,0,0,5,4]
  so4 = symmetryOrbits(BigFloat,Triangle())[7]
  positions = BigFloat[]
  weights = BigFloat[]

  for i in axes(xy,1)
    λ = barycentricCoordinates(refIn, xy[i,:])
    push!(positions, so4.compact(λ)...)
    push!(weights, wxy[i]*2)
  end

  for i in eachindex(t)
    λ = barycentricCoordinates(refIn, SVector{2,BigFloat}((t[i],t[i])))
    push!(positions, λ...)
    push!(weights, wt[i]*2)
  end

  cqr = CompactQuadratureRule(Triangle(), 7, orbits, positions)
  qr = expand(cqr)
  println(length(qr.points))
  @assert testWeights(qr; tol=sqrt(eps(Float64)))
  @assert testQuadratureRule(qr; tol=sqrt(eps(Float64)))

  ocqr = optimize(cqr; show_trace=true)
  oqr = expand(ocqr)
  @assert testWeights(oqr)
  @assert testQuadratureRule(oqr)

  write_file(file, ocqr; reference="Gun75", precision=32)
end
