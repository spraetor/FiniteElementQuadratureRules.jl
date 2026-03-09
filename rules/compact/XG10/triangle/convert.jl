using FiniteElementQuadratureRules
using FiniteElementQuadratureRules: barycentricCoordinates

function transformCompactRule3(file)
    x=[ 0.1551275866992061,
        0.9252314250365582,
        0.1384253518864596,
        0.3609812845058779,
        0.5812005281770837,
        0.3458850891448130,
        0.4510519110813880,
        0.6893778940829133,
        0.7634762995170715,
        0.2563805151012611,
        0.0,0.0,0.0,0.0]
    y=[ 0.3926587334507982,
        -0.5241973822167591,
        0.8231152771873217,
        0.0267765158868186,
        0.0484050181187623,
        0.4561437111102405,
        -0.3059662347573659,
        -0.5226600964300047,
        -0.3046970738749675,
        -0.5238003914950184,
        0.7059408419246426,
        1.0463216388186550,
        0.0240306520799389,
        -0.3100722765497620,]
    # w=[ 0.05567338455800167,
    #     0.00896962036448401,
    #     0.02240118142259882,
    #     0.06089917117558221,
    #     0.03108755935951332,
    #     0.03146650526166210,
    #     0.06246416621868210,
    #     0.02582709676771077,
    #     0.03178319217145733,
    #     0.03640547863919981,
    #     0.03498602106696763,
    #     0.01202119786551837,
    #     0.07686209708916766,
    #     0.07278083120266070,]
    ref1=ReferenceElement(Triangle, SVector{2,BigFloat}[[-1,-1/sqrt(big(3))], [0,2/sqrt(big(3))], [1,-1/sqrt(big(3))]],
      [[1,2], [1,3], [2,3]], 3/sqrt(big(3)))
    ref2=ReferenceElement(Triangle, SVector{2,BigFloat}[[-1,-1], [1,-1], [-1,1]],
      [[1,2], [1,3], [2,3]], BigFloat(2))
    geo=MultiLinearGeometry(ref1,ref2.coordinates)
    coordinates=map(xy->barycentricCoordinates(Triangle(), geo(big.([xy[1],xy[2]]))), zip(x,y))
    orbits=[0,0,0,10,0,4]
    positions=BigFloat[]
    so = symmetryOrbits(BigFloat, Triangle())
    j = 1
    for i in eachindex(orbits)
      for _ in 1:orbits[i]
        push!(positions, so[i].compact(coordinates[j])...)
        j += 1
      end
    end
    cqr = CompactQuadratureRule(Triangle(), 10, orbits, positions)
    qr = expand(cqr)
    @assert testWeights(qr; tol=BigFloat(sqrt(eps(Float64))))
    @assert testQuadratureRule(qr; tol=BigFloat(sqrt(eps(Float64))))

    ocqr = optimize(cqr; show_trace=true)
    println(typeof(ocqr))
    oqr = expand(ocqr)
    println(typeof(oqr))
    @assert testWeights(oqr)
    @assert testQuadratureRule(oqr)

    weights = BigFloat[]
    j = 1
    for i in eachindex(orbits)
      for _ in 1:orbits[i]
        push!(weights, oqr.weights[j])
        j += so[i].size
      end
    end
    weights=map(wi->wi/2, weights)

    cqr2 = CompactQuadratureRuleWithWeights(Triangle(), 10, orbits, ocqr.positions, weights)
    qr2 = expand(cqr2)
    @assert testWeights(qr2)
    @assert testQuadratureRule(qr2)

    write_file(file, cqr2; reference="XG10", precision=32)
  end
