include("InputParameters.jl")
t_e = 1800
t_m = 1200
include("InputVariables.jl")
include("ModuleSetup.jl")

#=-----------------
Simulation 
------------------=#
t_span = (0, 3600.0);
p = [
    # atmospheric and soil parameters
    soil_fine.α => finesoil.α, soil_fine.n => finesoil.n,
    soil_fine.Kₛ => finesoil.K_s, soil_fine.l => finesoil.l,
    soil_fine.θₛ => finesoil.θ_s, soil_fine.θᵣ => finesoil.θ_r,
    soil_fine.dz => Δz,
    # root parameters
    wup.εₓ => εₓ, wup.rᵣ => r_root, wup.kᵣ => k_r, wup.kₓ => k_x
]
u0 = [
    [soil_fine.h[i] => h0[i] for i in 1:nz]...,    # initial pressure head [cm or m]
    [wup.Hₓ[i] => h0[1] - z[1] for i in 1:nz]...,
    [wup.uptake[i] => 0.0 for i in 1:nz]..., 
]

guess_H = [
    [wup.Hᵣₛ[i] => h0[1] - z[1] for i in 1:nz]...,
]

ode_prob = ODEProblem(system_base, merge(Dict(u0), Dict(p)), t_span; guesses = Dict(guess_H))

@time sol_ode_3 = solve(ode_prob, Rodas5(), reltol=1e-5, abstol=1e-6, saveat = 1.0)
