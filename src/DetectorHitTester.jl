struct DetectorHitTester
    scaleFactor::Float64
    emin::Float64
    emax::Float64
    implement_epotek_cutoff::Bool
    implement_imperfect_mirror::Bool
    implement_QE::Bool
    implement_CE::Bool
end

function DetectorHitTester(;
    scaleFactor::Float64=1.0,
    emin::Float64=1.75,
    emax::Float64=7.00,
    implement_epotek_cutoff::Bool=true,
    implement_imperfect_mirror::Bool=true,
    implement_QE::Bool=true,
    implement_CE::Bool=true,
)
    DetectorHitTester(
        scaleFactor,
        emin,
        emax,
        implement_epotek_cutoff,
        implement_imperfect_mirror,
        implement_QE,
        implement_CE,
    )
end

function efficiency(hitTester::DetectorHitTester, energy::Float64; glueLayers::Int=1)::Float64
    if energy < dht.emin || energy > dht.emax
        return 0.0
    end
    result = 1.0
    if hitTester.implement_QE
        result *= QE_interpolator(energy)
    end
    if hitTester.implement_CE
        result *= CE()
    end
    if hitTester.implement_epotek_cutoff
        result *= epotek_305_interpolator(energy)^glueLayers
    end
    if hitTester.implement_imperfect_mirror
        result *= mirror_reflect(energy)
    end
    return result
end

function testHit(hitTester::DetectorHitTester, energy::Float64; glueLayers::Int=1)::Bool
    random_val = rand() * hitTester.scaleFactor
    return random_val < efficiency(hitTester, energy, glueLayers=glueLayers)
end
