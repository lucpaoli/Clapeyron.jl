# function VT_volume_function_coefficient(model, V, T, z)
#     μ_res = VT_chemical_potential_res(model, V, T, z)
#     Φ = @. 1 / exp(μ_res / R̄ / T)
#     return Φ
# end

# Mikyska V,T stability analysis
function VT_chemical_stability_analysis(model, V, T, z)
    cz = z ./ V
    μ(x) = VT_chemical_potential(model, 1.0, T, x)
    p(x) = pressure(model, 1.0, T, x)

    # Tangent lane distance function in concentration space
    tdp(cx) = sum((μ(cx) .- μ(cz) .* cx .- (p(cx) .- p(cz))))

    pure_models = split_model(model)
    psat = map(i -> i[1], saturation_pressure.(pure_models, T))

    # Trial vapourlike phase
    p0_vap = sum(psat .* z)
    x0_vap = psat .* z ./ p0_vap
    V0_vap = volume(model, p0_vap, T, x0_vap)
    cx0_vap = x0_vap ./ V0_vap

    # Trial liquidlike phase
    x0_liq = z ./ psat ./ (sum(z ./ psat))
    p0_liq = sum(psat .* x0_liq)
    V0_liq = volume(model, p0_liq, T, x0_liq)
    cx0_liq = x0_liq ./ V0_liq

    # function F_obj!(Fx, x)
    #     Fx = tdp(x .^ 2)
    #     return Fx
    # end
    # function J_obj!(Jx, x)
    #     Jx = ForwardDiff.jacobian(F_obj!, x)
    #     return Jx
    # end
    # function FJ_obj!(Fx, Jx, x)
    #     F_obj!(Fx, x)
    #     J_obj!(Jx, x)
    #     return Fx, Jx
    # end
    # function Jvop_obj!(x)
    #     function JacV(Fv, v)
    #             ForwardDiff.jacobian(JCache,F_obj!,)
    #     end
    #     return LinearMap(JacV, length(x))
    # end
    # vectorobj = NLSolvers.VectorObjective(F_obj!, J_obj!, FJ_obj!, Jvop_obj!)



    # tdp_func(cx0) = Solvers.optimize(cx -> tdp(cx .^ 2), sqrt.(cx0))
    tdp_func(cx0) = Solvers.optimize(cx -> tdp(cx), cx0)
    # res_vec = [tdp_func(cx0_liq), tdp_func(cx0_vap)]
    res_vec = tdp_func(cx0_vap)
    return res_vec
    # (tdp_min, idx) = findmin(map(x -> x.info.minimum, res_vec))
    # tdp_xmin = map(x -> sqrt.(abs.(x.info.solution)), res_vec)[idx]

    # return (tdp_xmin, tdp_min)
    # Φ(x) = VT_volume_function_coefficient(model, V, T, x)
    # Rewrite stability conditions in log space, using volume fractions
    # Rᴺ → Rᴺ
    # tdp_stationary_obj(cx) = log.(cx ./ cz) - log.(Φ(cz)) - log.(Φ(cx))

    # function f!(f, cx)
    #     f = tdp_stationary_obj(cx .^ 2)
    # end
    # tdp_func(cx0) = Solvers.nlsolve(f!, cx0)

    # res_vec = [tdp_func(cx0_liq), tdp_func(cx0_vap)]
    # tdp_xmin = map(x -> sqrt.(abs.(x.info.solution)), res_vec)
    # tdp_min = map(x -> tdp(x), tdp_xmin)
    # return (tdp_xmin, tdp_min)
end
