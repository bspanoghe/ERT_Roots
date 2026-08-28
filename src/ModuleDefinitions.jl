@independent_variables t
D = Differential(t)

# # Modules
# ## Translated
function soil_module(; name, Ψ_m, α, n, Kₛ, l, θₛ, θᵣ, dz, z)
    ρ_w = 1.0 # g / cm^3
    g = 9.8 * 1.0e-5 # hN / g
    Pₕ = ρ_w * g * z # MPa

    params = @parameters(
        α = α, [description = ""],
        n = n, [description = ""],
        Kₛ = Kₛ, [description = ""],
        l = l, [description = ""],
        θₛ = θₛ, [description = ""],
        θᵣ = θᵣ, [description = ""],
        dz = dz, [description = "Layer width"],
        Pₕ = Pₕ, [description = "Gravitational water potential"],
    )
    vars = @variables (
        Ψ(t), [description = "Total water potential [?]"], # alias `hT`
        Ψ_m(t) = Ψ_m, [description = "Matrix water potential [?]"], # alias `h`
        C(t), [description = "Soil water capacitance [?]"],
        K(t), [description = "Hydraulic conductivity [?]"],
        θ(t), [description = "Volumetric water content [?]"],
        # s(t), [description = "Root water uptake sink term (?) [?]"],
        F(t), [description = "Water flux [?]"], # alias `q`
        ΣF(t), [description = "Net water influx [?]"], # alias `dq`
    )
    eqs = [
        C ~ vanGenuchten_C(Ψ, θₛ, θᵣ, α, n),
        K ~ vanGenuchten_K(Ψ, θₛ, θᵣ, α, n, Kₛ, l),
        θ ~ vanGenuchten_θ(Ψ, θₛ, θᵣ, α, n),

        D(Ψ_m) ~ ( ΣF/dz #= - s/dz =# ) / C,
        Ψ ~ Ψ_m + Pₕ # eq. 7 in paper
    ]

    system = ODESystem(eqs, t; name)
    return system
end

# ## TODO

function rootuptake_module(; name,  εₓ, rᵣ, kᵣ, kₓ, Ψ_ref, k_Ψ, kc)
    params = @parameters (εₓ = εₓ, rᵣ = rᵣ, kᵣ = kᵣ, kₓ = kₓ, Ψ_ref = Ψ_ref, k_Ψ = k_Ψ, kc = kc)
    vars = @variables (        
        F(t)[1:nz], [description = ""],
        Tp(t), [description = ""],
        Kₓ(t)[1:nz], [description = ""],
        dz(t), [description = ""],
        lᵣ(t)[1:nz], [description = ""],
        Hₓ(t)[1:nz], [description = ""],
        dHₓ(t)[1:nz], [description = ""],
        dWₓ(t)[1:nz], [description = ""],
        Wₓ(t)[1:nz], [description = ""],
        hₛ(t)[1:nz], [description = ""],
        dhₛ(t), [description = ""],
        Hₛ(t)[1:nz], [description = ""],
        s(t)[1:nz], 
        r_rhiz(t)[1:nz], [description = ""],
        rld(t)[1:nz], 
        ρ(t)[1:nz], [description = ""],
        B(t)[1:nz], [description = ""],
        Hᵣₛ(t)[1:nz], [description = ""],
        Kᵣ(t)[1:nz], 
        uptake(t)[1:nz], 
        H₀(t), [description = ""],
        LAI(t), [description = ""],
        f_Ψ(t), [description = ""],
    )
    eqs = [
        # root-soil dimension aspects
        [lᵣ[i] ~ rld[i] * dz for i in 1:nz]..., # implicit multiplication with 1 cm² soil surface area
        [r_rhiz[i] ~ rhizo_fun(rld[i]) for i in 1:nz]...,
        [ρ[i] ~ r_rhiz[i] / rᵣ for i in 1:nz]...,
        [B[i] ~ B_fun(ρ[i]) for i in 1:nz]...,

        # Hydraulic conductivities based on root dimensions and intrinsic k's (h⁻¹)
        [Kᵣ[i] ~  kᵣ * (2*π * rᵣ * lᵣ[i]) / rᵣ  for i in 1:nz]...,
        [Kₓ[i] ~ kₓ / dz * (lᵣ[i]/dz * π * rᵣ^2) for i in 1:nz]...,

        # root xylem flows (positive upward) (cm h⁻¹)
        F[1] ~ (Hₓ[1] - H₀)*Kₓ[1],
        [F[i] ~ (Hₓ[i] - Hₓ[i-1])*(Kₓ[i]*Kₓ[i-1]*2/(Kₓ[i] + Kₓ[i-1])) for i in 2:nz]..., 
        
        # solving for root-soil interface water potential (cm)
        [kᵣ * Hₓ[i] + B[i] * Hₛ[i] * ksoilfun(hₛ[i], Hᵣₛ[i] - z[i]) ~
         Hᵣₛ[i] * (B[i] * ksoilfun(hₛ[i], Hᵣₛ[i] - z[i]) + kᵣ)
         for i in 1:nz]...,

        # water uptake s (cm h⁻¹) per soil layer i
        [s[i] ~ (Hᵣₛ[i] - Hₓ[i]) * Kᵣ[i] for i in 1:nz]...,

        # representative "water mass" of the root xylem (g)
        [Wₓ[i] ~ lᵣ[i] * π * ρ_w * rᵣ^2 for i in 1: nz]...,

        # change in water mass of the root xylem (only for change in Hₓ, not for integration!!!)
        dWₓ[nz] ~ s[nz] - F[nz], 
        [dWₓ[i] ~ s[i] - F[i] + F[i+1] for i in 1:(nz-1)]...,

        # change in xylem water potential (cm h⁻¹)
        dHₓ[1] ~ dh_fun(Wₓ[1], dWₓ[1], εₓ, tiny, dhₛ),
        [dHₓ[i] ~ dh_fun(Wₓ[i], dWₓ[i], εₓ, tiny, dHₓ[i-1])  for i in 2:nz]..., 

        # collar water potential (cm)
        H₀ ~ max(Hₓ[1] - Tp/(Kₓ[1] + tiny), hmin),
        #Ψ ~ H₀ * 98.1e-6,
        f_Ψ ~ 0.5 *( 1 + tanh(k_Ψ * (H₀ * 98.1e-6 - Ψ_ref))),
        Tp ~ ET0_inputfun(t) /10.0 * kc * (1 - exp(-0.45*(LAI))) * f_Ψ,
        # differential equations
        [D(Hₓ[i]) ~ dHₓ[i] for i in 1:nz]...,
        [D(uptake[i]) ~ s[i] for i in 1:nz]...,
    ]
    system = ODESystem(eqs, t; name)
    return system
end

function phenology_module(; name, Tmin, Tmax, Topt, v_max, S_ref, k_s, k_Ψ, Ψ_ref, r_LAI, r_max)
    params = @parameters(
        Tmin = Tmin, [description = ""],
        Tmax = Tmax, [description = ""],
        Topt = Topt, [description = ""],
        v_max = v_max, [description = ""],
        S_ref = S_ref, [description = ""],
        k_s = k_s, [description = ""],
        k_Ψ = k_Ψ, [description = ""],
        Ψ_ref = Ψ_ref, [description = ""],
        r_LAI = r_LAI, [description = ""],
        r_max = r_max, [description = ""],
    )
    vars = @variables (
        f_T(t), [description = ""],
        f_R(t), [description = ""],
        f_Ψ(t), [description = ""],
        T(t), [description = ""],
        Sᵥ(t), [description = ""],
        dSᵥ(t), [description = ""],
        LAI(t), [description = ""],
        dLAI(t), [description = ""],
        Ψ(t), [description = ""],
        Sᵣ(t), [description = ""],
        dSᵣ(t), [description = ""],
    )
    eqs = [
        f_T ~ (((Tmax - T)/(Tmax-Topt))*((T - Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)))^1.0,
        f_R ~ 0.5 *( 1 + tanh(k_s * (S_ref - Sᵥ))),
        f_Ψ ~ 0.5 *( 1 + tanh(k_Ψ * (Ψ - Ψ_ref))),
        dSᵥ ~ v_max * f_T * f_R * 0.5 *( 1 + tanh(10 * (t - 1000))),
        dSᵣ ~ (1-f_R) * r_max * f_T ,
        dLAI ~ dSᵥ * r_LAI * Sᵥ * (S_ref - Sᵥ)/S_ref * f_Ψ,
        D(Sᵥ) ~ dSᵥ,
        D(LAI) ~ dLAI,
        D(Sᵣ) ~ dSᵣ, 
    ]
    system = ODESystem(eqs, t; name)
    return system
end

# # Module connections
function soil_connection(; name, dz)
    params = @parameters(
        dz = dz, [description = "Layer width"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1"],
        K_half(t), [description = "Hydraulic conductivity of connection"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2"],
        Ψ_1(t), [description = "Total water potential of compartment 1"],
        Ψ_2(t), [description = "Total water potential of compartment 2"],
    )
    eqs = [
        F ~ K_half * ( (Ψ_2-Ψ_1) / dz -  dz/dz ),
        K_half ~ 2 / (1/K_1 + 1/K_2),
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = [
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
        connection_MTK.K_1 ~ node_MTK.K,
        connection_MTK.K_2 ~ nb_node_MTK.K,
    ]
    return System(eqs, t; name), get_connection_eqset
end