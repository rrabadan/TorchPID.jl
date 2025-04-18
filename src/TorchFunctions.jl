"""
`SellmeierCoefficients` is a struct representing Sellmeier coefficients for refractive index calculations.
    The Sellmeier equation is used to model the refractive index as a function of wavelength. 
    The equation has the form: n²(λ) = 1 + A₁λ²/(λ²-B₁) + A₂λ²/(λ²-B₂) + A₃λ²/(λ²-B₃), 
    Where our coefficients map as: a=A₁, b=B₁, c=A₂, d=B₂, e=A₃, f=B₃

# Fields
- `a::Float64`: First numerator coefficient in the Sellmeier equation.
- `b::Float64`: First denominator coefficient in the Sellmeier equation.
- `c::Float64`: Second numerator coefficient in the Sellmeier equation.
- `d::Float64`: Second denominator coefficient in the Sellmeier equation.
- `e::Float64`: Third numerator coefficient in the Sellmeier equation.
- `f::Float64`: Third denominator coefficient in the Sellmeier equation.
"""
struct SellmeierCoefficients
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    e::Float64
    f::Float64
end


"""
`CORNING` is a constant containing the Sellmeier coefficients for Corning 7980 fused silica (fused quartz).
These coefficients are used to calculate the refractive index as a function of wavelength.
"""
const CORNING = SellmeierCoefficients(
    0.68374049400,
    0.00460352869,
    0.42032361300,
    0.01339688560,
    0.58502748000,
    64.49327320000,
)

"""
    ngroup_Corning(energy_input::Float64)
    ngroup_Corning(energy_input::Float64, nphase::Float64)

`ngroup_Corning` calculates the group refractive index of Corning material (quartz). 
    The group refractive index quantifies how the group velocity of light varies within the material. 
    It is derived from the phase refractive index using the relation: nₗ = n - λ(dn/dλ). 
    This calculation employs the Sellmeier equation and its derivative, with the wavelength in μm computed from energy in eV using the formula: λ = 1.24/E.

# Arguments
- `energy_input::Float64`: Photon energy in eV.
- `nphase::Float64`: Phase refractive index at the given energy.

# Returns
- `Float64`: The group refractive index value.
"""
function ngroup_Corning(energy_input)
    lambda_um = 1.24 / energy_input
    lambda_um_sq = lambda_um^2

    n_sq =
        1 +
        (CORNING.a * lambda_um_sq) / (lambda_um_sq - CORNING.b) +
        (CORNING.c * lambda_um_sq) / (lambda_um_sq - CORNING.d) +
        (CORNING.e * lambda_um_sq) / (lambda_um_sq - CORNING.f)

    nphase = sqrt(n_sq)

    dndlambda =
        (lambda_um / nphase) * (
            (
                CORNING.a / (lambda_um_sq - CORNING.b) -
                (CORNING.a * lambda_um_sq) / ((lambda_um_sq - CORNING.b)^2)
            ) +
            (
                CORNING.c / (lambda_um_sq - CORNING.d) -
                (CORNING.c * lambda_um_sq) / ((lambda_um_sq - CORNING.d)^2)
            ) +
            (
                CORNING.e / (lambda_um_sq - CORNING.f) -
                (CORNING.e * lambda_um_sq) / ((lambda_um_sq - CORNING.f)^2)
            )
        )

    ngroup = nphase - lambda_um * dndlambda

    return ngroup
end

function ngroup_Corning(energy_input::Float64, nphase::Float64)::Float64
    lambda_um = 1.24 / energy_input
    lambda_um_sq = lambda_um^2

    dndlambda =
        (lambda_um / nphase) * (
            (
                CORNING.a / (lambda_um_sq - CORNING.b) -
                (CORNING.a * lambda_um_sq) / ((lambda_um_sq - CORNING.b)^2)
            ) +
            (
                CORNING.c / (lambda_um_sq - CORNING.d) -
                (CORNING.c * lambda_um_sq) / ((lambda_um_sq - CORNING.d)^2)
            ) +
            (
                CORNING.e / (lambda_um_sq - CORNING.f) -
                (CORNING.e * lambda_um_sq) / ((lambda_um_sq - CORNING.f)^2)
            )
        )

    ngroup = nphase - lambda_um * dndlambda

    return ngroup
end

"""
    nphase_Corning(energy_input::Float64)

`nphase_corning` calculates the phase refractive index of Corning material. 
    The phase refractive index quantifies how the phase velocity of light propagates through the material.
    It is calculated using the Sellmeier equation with the CORNING coefficients. 
    The wavelength in μm is derived from the photon energy in eV using the formula: λ = 1.24/E.

# Arguments
- `energy_input::Float64`: Photon energy in eV.

# Returns
- `Float64`: The phase refractive index value.
"""
function nphase_Corning(energy_input::Float64)::Float64
    lambda_um = 1.24 / energy_input
    lambda_um_sq = lambda_um^2

    n_sq =
        1 +
        (CORNING.a * lambda_um_sq) / (lambda_um_sq - CORNING.b) +
        (CORNING.c * lambda_um_sq) / (lambda_um_sq - CORNING.d) +
        (CORNING.e * lambda_um_sq) / (lambda_um_sq - CORNING.f)

    return sqrt(n_sq)
end

"""
   CE()

`CE` returns the collection efficiency constant for the detector system.
    This represents the fraction of Cherenkov photons that are successfully collected
    by the optical system. This value accounts for various geometric and optical factors in the detector design.

# Returns
- `Float64`: Collection efficiency value (0.65).
"""
function CE()
    return 0.65
end

"""
    epotek_305_interpolator(energy_input::Float64)

`epotek_305_interpolator` interpolates the transmission probability through Epotek-305 optical coupling material.
    Valid for wavelengths between 200nm and 800nm. Returns 0 for wavelengths outside this range.
    Wavelength is calculated from energy using the formula: λ(nm) = 1240/E(eV)

# Arguments
- `energy_input::Float64`: Photon energy in eV.

# Returns
- `Float64`: Transmission probability value between 0 and 1.
"""
function epotek_305_interpolator(energy_input::Float64)
    lambda = 1240.0 / energy_input
    if lambda < 200.0 || lambda > 800.0
        return 0.0
    end

    T_epotek_305 = [
        0.998881599,
        0.998828541,
        0.998902983,
        0.998924273,
        0.998934765,
        0.998902995,
        0.99888161,
        0.999030512,
        0.998881587,
        0.999009342,
        0.999094531,
        0.999115751,
        0.999019817,
        0.999041258,
        0.998870634,
        0.998838968,
        0.998679277,
        0.998934686,
        0.999073216,
        0.998828279,
        0.998923987,
        0.998881468,
        0.998870935,
        0.998711148,
        0.998828242,
        0.99893472,
        0.99890282,
        0.998966601,
        0.998923952,
        0.998860162,
        0.998849313,
        0.998977146,
        0.998966568,
        0.998902726,
        0.998913252,
        0.998977287,
        0.998796195,
        0.998859992,
        0.998817325,
        0.998860114,
        0.998923734,
        0.998934413,
        0.998753529,
        0.99899838,
        0.998849338,
        0.998838486,
        0.998891884,
        0.998934561,
        0.998881277,
        0.998849301,
        0.999009025,
        0.998966425,
        0.998934436,
        0.99883856,
        0.998934402,
        0.998870695,
        0.998838585,
        0.998987746,
        0.998732127,
        0.998753449,
        0.998859919,
        0.998891801,
        0.998859786,
        0.998881336,
        0.999083644,
        0.998902644,
        0.99896648,
        0.99896648,
        0.999179437,
        0.999392363,
        0.999168381,
        0.99918967,
        0.99911512,
        0.999147131,
        0.999040614,
        0.999029975,
        0.999030016,
        0.998912893,
        0.998934186,
        0.998891471,
        0.998891601,
        0.998934084,
        0.998817035,
        0.9989022,
        0.999029923,
        0.998944754,
        0.999051335,
        0.998902118,
        0.998934141,
        0.99885953,
        0.999051365,
        0.999104602,
        0.999072673,
        0.998998049,
        0.998944844,
        0.99890227,
        0.998838226,
        0.998902106,
        0.998955335,
        0.99892347,
        0.998880681,
        0.998912731,
        0.999051112,
        0.998944653,
        0.998997943,
        0.998923287,
        0.99895538,
        0.998933981,
        0.9988914,
        0.999104497,
        0.998859397,
        0.998987293,
        0.998933902,
        0.998933868,
        0.998912557,
        0.998997772,
        0.998923195,
        0.998997825,
        0.998997943,
        0.99881653,
        0.998795515,
        0.998859324,
        0.999008455,
        0.998933913,
        0.999051021,
        0.998901907,
        0.998901849,
        0.99882723,
        0.998901884,
        0.999061714,
        0.998955246,
        0.998859214,
        0.998880502,
        0.998976382,
        0.998944383,
        0.998859178,
        0.998912499,
        0.999061614,
        0.999008201,
        0.99895499,
        0.998955068,
        0.999072238,
        0.998944372,
        0.998923057,
        0.998891045,
        0.998986904,
        0.999104353,
        0.998976317,
        0.998976262,
        0.998891128,
        0.99891222,
        0.99882718,
        0.998912406,
        0.999029478,
        0.99896562,
        0.998987012,
        0.998976295,
        0.998986861,
        0.999040276,
        0.998869623,
        0.998933686,
        0.998741843,
        0.998880311,
        0.998997536,
        0.999040102,
        0.998922793,
        0.998933493,
        0.998922782,
        0.998912174,
        0.9989654,
        0.999029302,
        0.998880179,
        0.999114789,
        0.998944023,
        0.998997333,
        0.999018604,
        0.99904002,
        0.998848184,
        0.998944101,
        0.999093401,
        0.998858715,
        0.998965278,
        0.999103885,
        0.998880048,
        0.998869273,
        0.998944124,
        0.999189307,
        0.999103885,
        0.998976,
        0.998943966,
        0.999061133,
        0.998954656,
        0.998911942,
        0.998943865,
        0.999007746,
        0.998911907,
        0.998879845,
        0.998943763,
        0.998826504,
        0.999061153,
        0.998986483,
        0.998954332,
        0.998986472,
        0.998922506,
        0.998965079,
        0.998890418,
        0.998965146,
        0.99893322,
        0.998965068,
        0.998943673,
        0.998901204,
        0.998975618,
        0.998826416,
        0.999007725,
        0.998858399,
        0.998954466,
        0.998922402,
        0.998965024,
        0.998932992,
        0.998879665,
        0.998943527,
        0.998772941,
        0.998943775,
        0.998911524,
        0.998911698,
        0.99887963,
        0.998794358,
        0.998890169,
        0.998964848,
        0.998836924,
        0.99878352,
        0.998783728,
        0.998804989,
        0.998975497,
        0.998911524,
        0.99905018,
        0.998943504,
        0.998932912,
        0.998900888,
        0.998943515,
        0.998868827,
        0.998922126,
        0.99886867,
        0.998762112,
        0.998794178,
        0.998932776,
        0.99884741,
        0.998890039,
        0.99888998,
        0.998900595,
        0.998825965,
        0.99887939,
        0.998836738,
        0.998911292,
        0.998804594,
        0.998857911,
        0.998964571,
        0.998836465,
        0.99882584,
        0.99881509,
        0.998836614,
        0.998708658,
        0.998900548,
        0.998953953,
        0.998900618,
        0.998911199,
        0.998825777,
        0.998847004,
        0.998921827,
        0.998740379,
        0.998932411,
        0.998654967,
        0.998708259,
        0.998793702,
        0.998697722,
        0.998783065,
        0.998857826,
        0.998900442,
        0.998814926,
        0.998783,
        0.998665613,
        0.998708259,
        0.998772286,
        0.998654867,
        0.998697514,
        0.998718825,
        0.998793624,
        0.998654895,
        0.998750814,
        0.998654838,
        0.998686753,
        0.998633326,
        0.998697458,
        0.998750747,
        0.998729447,
        0.998622884,
        0.998761451,
        0.998761226,
        0.998697235,
        0.998537231,
        0.998558524,
        0.998515956,
        0.998494635,
        0.998547986,
        0.998558678,
        0.998579712,
        0.998558432,
        0.998622648,
        0.99850512,
        0.998654608,
        0.998633311,
        0.998462639,
        0.998441253,
        0.998462327,
        0.998515718,
        0.998387682,
        0.998537106,
        0.998387631,
        0.998451516,
        0.998355614,
        0.998451368,
        0.998494088,
        0.99844092,
        0.998366031,
        0.998397932,
        0.998462114,
        0.998248799,
        0.998312597,
        0.998291365,
        0.998376606,
        0.998280355,
        0.998397932,
        0.998184459,
        0.998280355,
        0.998333974,
        0.998290982,
        0.998237669,
        0.998237631,
        0.998184168,
        0.99824805,
        0.998248443,
        0.998216269,
        0.99825877,
        0.99818409,
        0.998237669,
        0.998130781,
        0.998077211,
        0.998076985,
        0.998087729,
        0.998226534,
        0.998055722,
        0.998087954,
        0.998023547,
        0.998055618,
        0.99798105,
        0.997970151,
        0.997863636,
        0.998087566,
        0.997906318,
        0.997938012,
        0.99805541,
        0.997980424,
        0.997937968,
        0.99794863,
        0.997906095,
        0.997841811,
        0.997916377,
        0.997959205,
        0.997766665,
        0.99793777,
        0.997852472,
        0.997649548,
        0.997884322,
        0.997852289,
        0.99772414,
        0.997574735,
        0.997755835,
        0.997905535,
        0.997553314,
        0.997639192,
        0.997649322,
        0.997681475,
        0.997531656,
        0.997435432,
        0.997425269,
        0.997403707,
        0.997318663,
        0.99761747,
        0.997243973,
        0.99721136,
        0.997105659,
        0.997083832,
        0.997232965,
        0.997159152,
        0.997212194,
        0.997020249,
        0.99711597,
        0.997137456,
        0.997114706,
        0.997168077,
        0.997168047,
        0.997370671,
        0.997508368,
        0.997593712,
        0.997753434,
        0.997614284,
        0.997763677,
        0.997998651,
        0.997688678,
        0.997816337,
        0.997677381,
        0.997720487,
        0.997720462,
        0.99747472,
        0.997197022,
        0.997165442,
        0.997410629,
        0.997324801,
        0.997122379,
        0.996961787,
        0.996897699,
        0.996939769,
        0.996641999,
        0.996663137,
        0.996492429,
        0.996683995,
        0.996609771,
        0.996801763,
        0.996439607,
        0.996363753,
        0.995873776,
        0.995990034,
        0.996183247,
        0.995927181,
        0.995573802,
        0.995541061,
        0.995850134,
        0.995784924,
        0.995558505,
        0.995760534,
        0.99530039,
        0.99517234,
        0.994989937,
        0.994903531,
        0.994603508,
        0.995274932,
        0.994792169,
        0.994633216,
        0.994688654,
        0.994534172,
        0.994298514,
        0.99404092,
        0.994275239,
        0.994507852,
        0.993911589,
        0.994144835,
        0.993484155,
        0.993116315,
        0.993064488,
        0.992431471,
        0.990456298,
        0.989938819,
        0.990670323,
        0.990839155,
        0.989663889,
        0.989323042,
        0.989499426,
        0.989571723,
        0.988215488,
        0.988587854,
        0.988340663,
        0.988058581,
        0.987537436,
        0.986521725,
        0.986695564,
        0.986258019,
        0.984096706,
        0.984458735,
        0.985231547,
        0.985237202,
        0.9850423,
        0.985045807,
        0.984317215,
        0.984465019,
        0.984317977,
        0.983200214,
        0.982973066,
        0.983622206,
        0.982969962,
        0.982657286,
        0.982044599,
        0.982416993,
        0.982695322,
        0.982481564,
        0.981761222,
        0.981799451,
        0.981482674,
        0.981686437,
        0.981675702,
        0.980195153,
        0.981738738,
        0.980693489,
        0.980342926,
        0.980703843,
        0.980908945,
        0.981672443,
        0.979949794,
        0.980497653,
        0.980139911,
        0.980373641,
        0.980417194,
        0.980490632,
        0.980433406,
        0.979839921,
        0.978673647,
        0.979011975,
        0.979218495,
        0.978891084,
        0.978092949,
        0.978521516,
        0.978365488,
        0.977895867,
        0.977430369,
        0.976970999,
        0.976969762,
        0.976913078,
        0.976335747,
        0.976138805,
        0.975671871,
        0.974875074,
        0.975245077,
        0.974606077,
        0.974294165,
        0.972993321,
        0.972230886,
        0.971158504,
        0.970063489,
        0.968571398,
        0.966896047,
        0.965155743,
        0.962965757,
        0.959661258,
        0.955974436,
        0.952048814,
        0.947245233,
        0.941463994,
        0.934591643,
        0.926238271,
        0.916243819,
        0.90498153,
        0.891726008,
        0.876471415,
        0.858650405,
        0.838549214,
        0.815399429,
        0.789685031,
        0.760417456,
        0.727601438,
        0.691790099,
        0.652270042,
        0.609536921,
        0.564048063,
        0.515891784,
        0.466357006,
        0.415440737,
        0.364216905,
        0.314434726,
        0.267122394,
        0.222653363,
        0.181987253,
        0.145892721,
        0.114646009,
        0.08846426,
        0.067323606,
        0.050664488,
        0.037834717,
        0.028359364,
        0.021581714,
        0.016918647,
        0.013707002,
        0.011586385,
        0.010127606,
        0.009301563,
        0.008736251,
        0.008403361,
        0.008112813,
        0.00795355,
        0.007774038,
        0.007671203,
        0.007579419,
        0.007509125,
        0.007397785,
        0.007384792,
        0.007274694,
        0.007148539,
        0.007126338,
        0.007179635,
        0.007165711,
        0.007098442,
        0.007011981,
        0.007010913,
    ]

    index = floor(Int, 800 - lambda) + 1

    if index < 1 || index + 1 > 601
        return 0.0
    end

    e_left = 1240.0 / (800 - (index - 1))
    e_right = 1240.0 / (800 - index)
    delta = (energy_input - e_left) / (e_right - e_left)

    return T_epotek_305[index] + (T_epotek_305[index+1] - T_epotek_305[index]) * delta
end

"""
    mirrior_reflect(input_energy::Float64)

`mirror_reflect` calculates the mirror reflection probability for a given photon energy. 
    Applicable for wavelengths between 200nm and 800nm. Returns 0 for wavelengths outside this range. 
    The wavelength is derived from energy using the formula: λ(nm) = 1240/E(eV).
    The reflection probability is influenced by the mirror material and the angle of incidence.

# Arguments
- `input_energy::Float64`: Photon energy in eV.

# Returns
- `Float64`: Reflection probability value between 0 and 1.
"""
function mirror_reflect(input_energy::Float64)
    lambda = 1240.0 / input_energy
    if lambda < 200.0 || lambda > 800.0
        return 0.0
    end

    R_lambda = [
        0.860,
        0.870,
        0.876,
        0.880,
        0.8835,
        0.887,
        0.8895,
        0.891,
        0.8925,
        0.894,
        0.896,
        0.897,
        0.8975,
        0.8975,
        0.898,
        0.898,
        0.898,
        0.898,
        0.898,
        0.8975,
        0.897,
        0.8965,
        0.896,
        0.895,
        0.8945,
        0.8935,
        0.895,
        0.8945,
        0.893,
        0.8905,
        0.8915,
        0.892,
        0.8915,
        0.89,
        0.889,
        0.8875,
        0.886,
        0.8845,
        0.8835,
        0.882,
        0.88,
        0.8785,
        0.877,
        0.875,
        0.8735,
        0.8715,
        0.8695,
        0.867,
        0.8645,
        0.8625,
        0.8595,
        0.8565,
        0.854,
        0.85,
        0.8465,
        0.8425,
        0.838,
        0.833,
        0.8275,
        0.8225,
        0.817,
    ]

    index = floor(Int, (lambda - 200.0) / 10.0) + 1
    if index < 1 || index + 1 > length(R_lambda)
        return 0.0
    end

    e_upp = 1240.0 / (200.0 + 10.0 * (index - 1))
    e_low = 1240.0 / (210.0 + 10.0 * (index - 1))
    delta = (input_energy - e_low) / (e_upp - e_low)

    return R_lambda[index+1] + (R_lambda[index] - R_lambda[index+1]) * delta

end

"""
    QE_interpolator(input_energy::Float64)

`QE_interpolator` interpolates the quantum efficiency of the photon detector for a given photon energy. 
    Valid for wavelengths between 200nm and 800nm. Returns 0 for wavelengths outside this range.
    Wavelength is calculated from energy using the formula: λ(nm) = 1240/E(eV). 
    Quantum efficiency represents the probability that an incident photon generates a detectable electron in the photosensor.

# Arguments
- `input_energy::Float64`: Photon energy in eV.

# Returns
- `Float64`: Quantum efficiency value between 0 and 1.
"""
function QE_interpolator(input_energy::Float64)
    lambda = 1240.0 / input_energy
    if lambda < 200.0 || lambda > 800.0
        return 0.0
    end

    QE = [
        0.208,
        0.282,
        0.281,
        0.270,
        0.245,
        0.251,
        0.246,
        0.229,
        0.221,
        0.208,
        0.194,
        0.179,
        0.160,
        0.143,
        0.129,
        0.114,
        0.099,
        0.086,
        0.068,
        0.058,
        0.050,
        0.042,
        0.035,
        0.030,
        0.025,
        0.021,
        0.017,
        0.014,
        0.011,
        0.008,
        0.006,
    ]

    index = floor(Int, (lambda - 200.0) / 20.0) + 1

    if index < 1 || index + 1 > length(QE)
        return 0.0
    end

    e_upp = 1240.0 / (200.0 + 20.0 * (index - 1))
    e_low = 1240.0 / (220.0 + 20.0 * (index - 1))
    delta = (input_energy - e_low) / (e_upp - e_low)

    return QE[index+1] + (QE[index] - QE[index+1]) * delta
end
