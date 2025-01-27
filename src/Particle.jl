struct Particle
    pid::Int
    eventId::UInt
    trackId::UInt
    moduleId::UInt
    pMag::Float64
    pMagTrue::Float64
    xDir::Float64
    yDir::Float64
    zDir::Float64
    truePX::Float64
    truePY::Float64
    truePZ::Float64
    recoPX::Float64
    recoPY::Float64
    recoPZ::Float64
    recoPathlength::Float64
    truePathlength::Float64
    rotMatrix::Array{Float64,2}

    # Inner constructor with default values
    function Particle(; pid::Int=211, eventId::UInt=0, trackId::UInt=0, moduleId::UInt=0,
        pMag::Float64=0.0, pMagTrue::Float64=0.0,
        xDir::Float64=0.0, yDir::Float64=0.0, zDir::Float64=0.0,
        truePX::Float64=0.0, truePY::Float64=0.0, truePZ::Float64=0.0,
        recoPX::Float64=0.0, recoPY::Float64=0.0, recoPZ::Float64=0.0,
        recoPathlength::Float64=0.0, truePathlength::Float64=0.0,
        rotMatrix::Array{Float64,2}=Array{Float64,2}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]))
        new(pid, eventId, trackId, moduleId, pMag, pMagTrue, xDir, yDir, zDir, truePX, truePY, truePZ, recoPX, recoPY, recoPZ, recoPathlength, truePathlength, rotMatrix)
    end
end


function beta(p::Particle, m::Float64)
    return p.pMag / sqrt(p.pMag^2 + m^2)
end

function gamma(p::Particle, m::Float64)
    return 1 / sqrt(1 - beta(p, m)^2)
end

function initRotation(p::Particle)
    uperp::Float64 = sqrt(p.xDir^2 + p.yDir^2)

    if uperp != 0
        p.rotMatrix[1, 1] = p.xDir * p.zDir / uperp
        p.rotMatrix[1, 2] = -p.yDir / uperp
        p.rotMatrix[1, 3] = p.xDir
        p.rotMatrix[2, 1] = p.yDir * p.zDir / uperp
        p.rotMatrix[2, 2] = p.xDir / uperp
        p.rotMatrix[2, 3] = p.yDir
        p.rotMatrix[3, 1] = -uperp
        p.rotMatrix[3, 2] = 0
        p.rotMatrix[3, 3] = p.zDir
    end
end