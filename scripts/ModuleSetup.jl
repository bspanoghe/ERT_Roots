#=-----------------------------
Components of the model
-----------------------------=#

components = Vector{ODESystem}(undef, 0)

@independent_variables t, [description = "time (h)"]
D = Differential(t);
@constants (ρ_w = 1.0)

@named soil_fine = soil_module(α = soilpar.α, n = soilpar.n, Kₛ = soilpar.K_s, l = soilpar.l, θₛ = soilpar.θ_s, θᵣ = soilpar.θ_r,
        dz = Δz) 
@named wup = rootuptake_module(εₓ = εₓ, rᵣ = r_root, kᵣ = k_r, kₓ = k_x)


#=--------------------------------
Connections between the modules 
--------------------------------=#

connections = Vector{Expr}(undef, 0)

soil_fine_connections = :([
    soil_fine.prec_evap ~ - Precip_inputfun(t) + Evaporation(t)
])

wup_connections = :([
    wup.Tp ~ Transpiration(t)
    [wup.s[i] ~ soil_fine.s[i] for i in 1:nz]...
    [wup.Hₛ[i] ~ soil_fine.hT[i] for i in 1:nz]...
    [wup.hₛ[i] ~ soil_fine.h[i] for i in 1:nz]...
    [wup.rld[i] ~ root_fun(-z[i],t,t_m, t_e, d_e, d_s, rld_d ,p_max, p_half) for i in 1:nz]...
    wup.dz ~ soil_fine.dz
    wup.dhₛ ~ soil_fine.dh[1]
])


#=-----------------------------
Assemble the system 
------------------------------=#
push!(components, soil_fine)
push!(components, wup)
push!(connections, soil_fine_connections)
push!(connections, wup_connections)

evaluated_connections = reduce(vcat, eval.(connections))
system_base = ODESystem(evaluated_connections, t, name = :system, systems = components) |> structural_simplify 
