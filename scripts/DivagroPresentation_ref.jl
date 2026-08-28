include("InputParameters.jl")
include("InputRoots.jl")

sandyloam = (θ_s = 0.42, θ_r = 0.07, α  = 0.07, n  = 1.8, K_s = 15 /(24), l  = 0.5)
soilpar = sandyloam

include("InputSoil.jl")
##input

Ref_df = CSV.read("./data/meteo_ref.csv", DataFrame; delim = ',', header = true, ignorerepeated=true)
Ref_df.Timestamp = DateTime.(Ref_df.Timestamp, dateformat"yyyy-mm-ddTHH:MM:SS\Z")
Ref_df.Time = (Ref_df.Timestamp .- Ref_df.Timestamp[1]) ./ Hour(1)


efun = LinearInterpolation(Ref_df.Time, Ref_df.ETo_h)
ET0_inputfun(t) = efun(t) # mm h-¹
@register_symbolic ET0_inputfun(t)

pfun = linear_interpolation(Ref_df.Time, Ref_df.P)
Precip_inputfun(t) = pfun(t) # mm h-¹
@register_symbolic Precip_inputfun(t)
Precip_ref(t) = pfun(t)

tfun = linear_interpolation(Ref_df.Time, Ref_df.T)
T_inputfun(t) = tfun(t)
@register_symbolic T_inputfun(t)


components = Vector{ODESystem}(undef, 0)
@independent_variables t, [description = "time (h)"]
D = Differential(t);
@constants (ρ_w = 1.0)

@named soil_sl = soil_module(α = soilpar.α, n = soilpar.n, Kₛ = soilpar.K_s, l = soilpar.l, θₛ = soilpar.θ_s, θᵣ = soilpar.θ_r,
        dz = Δz) 
@named wup = rootuptake_module(εₓ = εₓ, rᵣ = r_root, kᵣ = k_r, kₓ = k_x, Ψ_ref = -1.0, k_Ψ = 5.0, kc = kc)
@named phenology = phenology_module(Tmin = 5.0, Tmax = 42.0, Topt = 25.0, v_max = 0.012, S_ref = 11.0, 
k_s = 1.0, k_Ψ = 4.0, Ψ_ref = -0.75, r_LAI = 0.5, r_max = 0.005)

#=--------------------------------
Connections between the modules 
--------------------------------=#

connections = Vector{Expr}(undef, 0)

soil_sl_connections = :([
    soil_sl.prec_evap ~ - Precip_inputfun(t)/10 + ET0_inputfun(t) /10* 1.0 * (1 - (1 - exp(-0.45*(phenology.LAI))))
])

wup_connections = :([
    wup.LAI ~ phenology.LAI  # from mm h-1 to cm d-1
    [wup.s[i] ~ soil_sl.s[i] for i in 1:nz]...
    [wup.Hₛ[i] ~ soil_sl.hT[i] for i in 1:nz]...
    [wup.hₛ[i] ~ soil_sl.h[i] for i in 1:nz]...
    [wup.rld[i] ~ root_fun(-z[i],t-800,t_m, t_e, d_e, d_s, rld_d ,p_max, p_half) for i in 1:nz]...
    wup.dz ~ soil_sl.dz
    wup.dhₛ ~ soil_sl.dh[1]
])

phenology_connections = :([
    phenology.Ψ ~ wup.H₀* 98.1e-6
    phenology.T ~ T_inputfun(t)
])

#=-----------------------------
Assemble the system 
------------------------------=#
push!(components, soil_sl)
push!(components, wup)
push!(components, phenology)

push!(connections, soil_sl_connections)
push!(connections, wup_connections)
push!(connections, phenology_connections)

evaluated_connections = reduce(vcat, eval.(connections))
system_base = ODESystem(evaluated_connections, t, name = :system, systems = components) |> structural_simplify 


t_span = (800, 4200.0);
u0 = [
    [soil_sl.h[i] => h0[i] for i in 1:nz]...,    # initial pressure head [cm or m]
    [wup.Hₓ[i] => h0[1] - z[1] for i in 1:nz]...,
    [wup.uptake[i] => 0.0 for i in 1:nz]..., 
    phenology.Sᵥ => 0.0, phenology.LAI => 0.0, phenology.Sᵣ => 0.0
]

guess_H = [
    [wup.Hᵣₛ[i] => h0[1] - z[1] for i in 1:nz]...,

]
p = Dict()

ode_prob = ODEProblem(system_base, merge(Dict(u0), Dict(p)), t_span; guesses = Dict(guess_H))

@time sol_ref = solve(ode_prob, Rodas5(), reltol=1e-5, abstol=1e-6, saveat = 1.0) # FBDF()


sol_df = DataFrame(
        Time = sol_ref.t,
        Sᵥ = sol_ref[phenology.Sᵥ],
        LAI = sol_ref[phenology.LAI],
        Sᵣ = sol_ref[phenology.Sᵣ],
        f_T = sol_ref[phenology.f_T],
        f_R = sol_ref[phenology.f_R],
        f_Ψ = sol_ref[phenology.f_Ψ],
        Ψ = sol_ref[phenology.Ψ],
        Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref.t))

plot(sol_df.Timestamp, sol_df.Sᵥ, label = "Sᵥ", linewidth = 2)
plot!(sol_df.Timestamp, sol_df.Sᵣ, label = "Sᵣ", linewidth = 2)

plot(sol_df.Timestamp, sol_df.LAI, label = "LAI", linewidth = 2)
plot(sol_df.Timestamp, sol_df.f_Ψ, label = "f_Ψ", linewidth = 2)
plot(sol_df.Timestamp, sol_df.Ψ, label = "Ψ", linewidth = 2)


plot(sol_ref, idxs=soil_sl.h, xlabel="Time", ylabel="h(t)", legend=:none)
plot(sol_ref, idxs=phenology.Ψ, xlabel="Time", ylabel="Ψ(t)", legend=:none)

plot(sol_ref, idxs=soil_sl.prec_evap, xlabel="Time", ylabel="Sᵥ(t)", legend=:none)
plot(sol_ref, idxs=wup.Tp, xlabel="Time", ylabel="f_Ψ(t)", legend=:none)

plot(sol_ref, idxs=wup.f_Ψ, xlabel="Time", ylabel="f_Ψ(t)", legend=:none)

plot(sol_ref, idxs=wup.rld, xlabel="Time", ylabel="f_Ψ(t)", legend=:none)