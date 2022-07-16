abstract type HSvdWModel <: vdWModel end

struct HSvdWParam <: EoSParam
    a::PairParam{Float64}
    b::PairParam{Float64}
    Tc::SingleParam{Float64}
    Pc::SingleParam{Float64}
    Vc::SingleParam{Float64}
    Mw::SingleParam{Float64}
end

struct HSvdW{T<:IdealModel,α,c,M} <: HSvdWModel
    components::Array{String,1}
    icomponents::UnitRange{Int}
    alpha::α
    mixing::M
    translation::c
    params::HSvdWParam
    idealmodel::T
    references::Array{String,1}
end
@registermodel HSvdW

export HSvdW

"""
    HSvdW(components::Vector{String};
    idealmodel=BasicIdeal,
    alpha = ClausiusAlpha,
    mixing = vdW1fRule,
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[], 
    ideal_userlocations=String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    verbose=false)

## Input parameters
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`
- `vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `k`: Pair Parameter (`Float64`)

## Model Parameters
- `a`: Pair Parameter (`Float64`)
- `b`: Pair Parameter (`Float64`)
- `Tc`: Single Parameter (`Float64`) - Critical Temperature `[K]`
- `Pc`: Single Parameter (`Float64`) - Critical Pressure `[Pa]`
- `Vc`: Single Parameter (`Float64`) - Molar Volume `[m^3/mol]`
- `Mw`: Single Parameter (`Float64`) - Molecular Weight `[g/mol]`


## Input models
- `idealmodel`: Ideal Model
- `alpha`: Alpha model
- `mixing`: Mixing model
- `activity`: Activity Model, used in the creation of the mixing model.
- `translation`: Translation Model

## Description


## References

Carnahan, N. F.; Starling, K. E. Intermolecular Repulsions and the Equation of State for Fluids. AIChE J. 1972, 18 (6), 1184–1189. https://doi.org/10.1002/aic.690180615.

"""
HSvdW

function HSvdW(components::Vector{String}; idealmodel=BasicIdeal,
    alpha=NoAlpha, # no alpha
    mixing=vdW1fRule, # 
    activity=nothing,
    translation=NoTranslation,
    userlocations=String[],
    ideal_userlocations=String[],
    alpha_userlocations=String[],
    mixing_userlocations=String[],
    activity_userlocations=String[],
    translation_userlocations=String[],
    verbose=false)
    params = getparams(components, ["properties/critical.csv", "properties/molarmass.csv", "SAFT/PCSAFT/PCSAFT_unlike.csv"]; userlocations=userlocations, verbose=verbose)
    k = params["k"]
    pc = params["pc"]
    Mw = params["Mw"]
    Tc = params["Tc"]
    Vc = params["vc"]
    init_mixing = init_model(mixing, components, activity, mixing_userlocations, activity_userlocations, verbose)
    a, b = ab_premixing(HSvdW, init_mixing, Tc, pc, Vc, k)
    init_idealmodel = init_model(idealmodel, components, ideal_userlocations, verbose)
    init_translation = init_model(translation, components, translation_userlocations, verbose)
    init_alpha = init_model(alpha, components, alpha_userlocations, verbose)
    icomponents = 1:length(components)
    packagedparams = HSvdWParam(a, b, Tc, pc, Vc, Mw)
    references = String["10.1002/aic.690180615"]
    model = HSvdW(components, icomponents, init_alpha, init_mixing, init_translation, packagedparams, init_idealmodel, references)
    return model
end

function ab_premixing(model::Type{<:HSvdWModel}, mixing::MixingRule, Tc, pc, vc, kij)
    _Tc = Tc.values
    # _Vc = vc.values
    _pc = pc.values
    components = vc.components
    a = epsilon_LorentzBerthelot(SingleParam("a", components, @. 0.359 * 1.38 * (R̄ * _Tc)^2 / _pc), kij)
    b = sigma_LorentzBerthelot(SingleParam("b", components, @. 0.359 * 0.5216 * R̄ * _Tc / _pc))
    return a, b
end

function ab_consts(::Type{<:HSvdWModel})
    Ωa = 0.359 * 1.38
    Ωb = 0.359 * 0.5216
    return Ωa, Ωb
end

# v = Z*RT/p
function a_res(model::HSvdWModel, V, T, z, _data=data(model, V, T, z))
    _, a, b, _ = _data

    RT = R̄ * T
    n = sum(z)

    v = V / n
    η = b / 4v
    Z_hs = (1 + η + η^2 - η^3) / (1 - η)^3 # "hard shell" Z

    # p = pressure(model, V, T, z)
    Z = (RT * v * (η^3 - η^2 - η - 1) + a * (η - 1)^3) / (RT * v * (η - 1)^3)

    return log(Z_hs / Z) - a / (R̄ * T * v)
end

