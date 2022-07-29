# One of the main issues with this (I think) is that it's not actually a cubic. It's defined as a derivative of vdW but is actually quintic
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

    Z = (RT * v * (η^3 - η^2 - η - 1) + a * (η - 1)^3) / (RT * v * (η - 1)^3)
    Z = Z_hs + a / (RT * v)

    # This is in molar quantities, should be corrected for n moles?
    return log(Z_hs / Z) - a * v / (RT)
end


function volume_impl(model::HSvdWModel, p, T, z=SA[1.0], phase=:unknown, threaded=true, vol0=nothing)
    #Threaded version
    TYPE = typeof(p + T + first(z))

    if ~isnothing(vol0)
        V0 = vol0
        V = _volume_compress(model, p, T, z, V0)
        return V
    end

    if phase != :unknown
        V0 = x0_volume(model, p, T, z, phase=phase)
        V = _volume_compress(model, p, T, z, V0)
        #isstable(model,V,T,z,phase) the user just wants that phase
        return V
    end
    if threaded
        Vg0 = x0_volume(model, p, T, z, phase=:v)
        Vl0 = x0_volume(model, p, T, z, phase=:l)
        _Vg = Threads.@spawn _volume_compress(model, $p, $T, $z, $Vg0)
        _Vl = Threads.@spawn _volume_compress(model, $p, $T, $z, $Vl0)
        Vg::TYPE = fetch(_Vg)
        Vl::TYPE = fetch(_Vl)
    else
        Vg0 = x0_volume(model, p, T, z, phase=:v)
        Vl0 = x0_volume(model, p, T, z, phase=:l)

        Vg = _volume_compress(model, p, T, z, Vg0)
        Vl = _volume_compress(model, p, T, z, Vl0)
    end

    #this catches the supercritical phase as well
    if isnan(Vl)
        _v_stable(model, Vg, T, z, phase)
        return Vg
    elseif isnan(Vg)
        _v_stable(model, Vl, T, z, phase)
        return Vl
    end

    err() = @error("model $model Failed to converge to a volume root at pressure p = $p [Pa], T = $T [K] and compositions = $z")
    if (isnan(Vl) & isnan(Vg))
        err()
        return zero(TYPE) / zero(TYPE)
    end
    _dfg, fg = ∂f(model, Vg, T, z)
    dVg, _ = _dfg
    gg = ifelse(abs((p + dVg) / p) > 0.03, zero(dVg) / one(dVg), fg + p * Vg)
    _dfl, fl = ∂f(model, Vl, T, z)
    dVl, _ = _dfl
    gl = ifelse(abs((p + dVl) / p) > 0.03, zero(dVl) / one(dVl), fl + p * Vl)
    if gg < gl
        _v_stable(model, Vg, T, z, phase)
        return Vg
    else
        _v_stable(model, Vl, T, z, phase)
        return Vl
    end
end
