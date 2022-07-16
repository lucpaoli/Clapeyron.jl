# Jindrova2015 - 10.1016/J.FLUID.2015.02.013
# This is a general Π-phase equilibrium calculation, meaning it does NOT assume a given number of phases a priori
# Operates within Helmholtz space
# Repeated VT-stability tests until a stable Π-phase state is found
# Uses Newton method with line-search, relying on modified Cholesky decompositon of Hessian. Guarantees decreasing values of total Helmholtz energy

# export tdp_minima

# function VT_volume_function_coefficient(model::EoSModel, V, T, z=SA[1.0])
#     μ_res = VT_chemical_potential_res(model, V, T, z)
#     Φ = @. 1 / exp(μ_res / R̄ / T)
#     return Φ
# end

# function tdp_stationary_obj!(F, model::EoSModel, V, T, cx, cz)
#     Φ(x) = VT_volume_function_coefficient(model, V, T, x)
#     # Rewrite stability conditions in log space, using volume fractions
#     # Rᴺ → Rᴺ
#     tdp_stationary_obj(cx) = @. log(cx / cz) - log(Φ(cz)) - log(Φ(cx))

#     F = tdp_stationary_obj(cx)
#     return F
# end

# # function tdp_stationary_grad!(∇F, model::EoSModel, V, T, cx, cz)
# #     Φ(x) = VT_volume_function_coefficient(model, V, T, x)
# #     tdp_stationary_obj(cx) = @. log(cx / cz) - log(Φ(cz)) - log(Φ(cx))

# #     ∇F = ForwardDiff.gradient(tdp_stationary_obj, cx)
# #     return ∇F
# # end

# # To guarantee location of all minima in the interval, IntervalRootFinding.jl may be useful
# function tdp_minima(model, V, T, z::AbstractVector{<:Number})
#     cz = z ./ V
#     μ(x) = VT_chemical_potential(model, 1.0, T, x)
#     p(x) = pressure(model, 1.0, T, x)

#     # Tangent lane distance function in concentration space
#     tdp(cx) = sum((μ(cx) .- μ(cz) .* cx .- (p(cx) .- p(cz))))

#     # ∂D/∂(cx) simplifies to μ' - μ
#     # Rᴺ → Rᴺ
#     # Does this have any solutions other than the trivial???
#     # tdp_stationary(cx) = μ(cx) .- μ(cz)

#     # Trial vapour phase
#     pure_models = split_model(model)
#     psat = [i[1] for i in saturation_pressure.(pure_models, T)]
#     p0 = sum(psat .* z)
#     x0 = psat .* z ./ p0
#     cx0 = x0 ./ V
#     # f!(F, x) = tdp_stationary_obj!(F, model, V, T, x, cz)
#     # ∇f!(∇F, x) = tdp_stationary_grad!(∇F, model, V, T, x, cz)
#     # fj!(F, ∇F, x) = (f!(F, x), ∇f!(∇F, x))
#     # fjop(F, ∇F, x) = (f!(F, x), ∇f!(∇F, x))

#     # vectorobj = NLSolvers.VectorObjective(f!, ∇f!,)
#     # stationary_points = nlsolve(f!, x0, autodiff=:forward)

#     Φ(x) = Clapeyron.VT_volume_function_coefficient(model, V, T, x)
#     # Rewrite stability conditions in log space, using volume fractions
#     # Rᴺ → Rᴺ
#     tdp_stationary_obj(cx) = log.(cx ./ cz) - log.(Φ(cz)) - log.(Φ(cx))
#     # stationary_points = Optim.optimize(tdp_stationary_obj, cx0, Optim.BFGS())
#     # return tdp_stationary_obj(cx0)
#     min_point = Optim.optimize(tdp, cx0, Optim.BFGS())
#     return min_point
# end

