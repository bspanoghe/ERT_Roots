@independent_variables t
D = Differential(t)

# # Modules
# ## Translated
function soil_module(; name, Ψ_m, α, n, Kₛ, l, θₛ, θᵣ, dz, z)
    ρ_w = 1.0, [description = "Density of water [g cm^-3]"], 
    g = 9.8 * 1.0e-5, [description = "Gravitational acceleration [MPa cm^2 g^-1]"],
    Pₕ = ρ_w * g * z, [description = "Gravitational water potential [MPa]"]

    params = @parameters(
        α = α, [description = "van Genuchten shape parameter (related to inverse of air entry suction, > 0) [cm^-1]"],
        n = n, [description = "van Genuchten shape parameter (related to pore size distribution, > 1) [-]"],
        Kₛ = Kₛ, [description = "Saturated hydraulic conductivity [cm h^-1]"],
        l = l, [description = "Pore connectivity parameter [-]"],
        θₛ = θₛ, [description = "Saturated volumetric water content [-]"],
        θᵣ = θᵣ, [description = "Residual volumetric water content [-]"],
        dz = dz, [description = "Layer width [cm]"],
        Pₕ = Pₕ, [description = "Gravitational water potential [MPa]"],
    )
    vars = @variables (
        Ψ(t), [description = "Total water potential [MPa]"], # alias `hT`
        Ψ_m(t), [description = "Matric water potential [MPa]"], # alias `h`
        h(t), [description = "Hydraulic head [cm]"], 
        C(t), [description = "Soil water capacitance [cm^-1]"],
        K(t), [description = "Hydraulic conductivity [cm h^-1]"], #eigenlijk moeten we dit zien als g per cm² per h, mits ρ_w = 1.0 g cm^-3
        θ(t), [description = "Volumetric water content [-]"],
        # s(t), [description = "Root water uptake sink term (?) [?]"],
        F(t), [description = "Water flux [cm h^-1]"], # alias `q` eigenlijk moeten we dit zien als g per cm² per h, mits ρ_w = 1.0 g cm^-3
        ΣF(t), [description = "Net water influx [cm h^-1]"], # alias `dq`  eigenlijk moeten we dit zien als g per cm² per h, mits ρ_w = 1.0 g cm^-3
    )
    eqs = [
        h * ρ_w * g ~ Ψ_m, 
        C ~ vanGenuchten_C(h, θₛ, θᵣ, α, n),
        K ~ vanGenuchten_K(h, θₛ, θᵣ, α, n, Kₛ, l),
        θ ~ vanGenuchten_θ(h, θₛ, θᵣ, α, n),

        D(h) ~ ( ΣF/dz #= - s/dz =# ) / C,
        Ψ ~ Ψ_m + Pₕ # eq. 7 in paper
    ]

    system = ODESystem(eqs, t; name)
    return system
end

# ## TODO

function rootuptake_module(; name,  εₓ, rᵣ, kᵣ, kₓ, Ψ_ref, k_Ψ, kc)
    params = @parameters (
        εₓ = εₓ, [description = "Root xylem elastic modulus [MPa]"],
        rᵣ = rᵣ, [description = "Root radius [cm]"],
        kᵣ = kᵣ, [description = "Intrinsic radial root hydraulic conductivity [h^-1]"],
        kₓ = kₓ, [description = "Intrinsic axial root hydraulic conductivity [h^-1]"],
        Ψ_ref = Ψ_ref, [description = "Reference water potential [MPa]"],
        k_Ψ = k_Ψ, [description = "Water potential sensitivity coefficient [-]"],
        kc = kc, [description = "Crop coefficient for transpiration [-]"]
        )
    vars = @variables (        
        F(t)[1:nz], [description = "Water flux from root xylem compartment i to root xylem compartment i-1 [cm h^-1]"],
        Tp(t), [description = "Transpiration rate [cm h^-1]"],
        Kₓ(t)[1:nz], [description = " Hydraulic conductivity of root xylem compartment i [cm h^-1]"],
        dz(t), [description = "Layer width [cm]"],
        lᵣ(t)[1:nz], [description = "Root length in soil layer i [cm]"],
        Hₓ(t)[1:nz], [description = "Water potential of root xylem compartment i [cm]"],
        dHₓ(t)[1:nz], [description = "Change in water potential of root xylem compartment i [cm h^-1]"],
        dWₓ(t)[1:nz], [description = "Change in water mass of root xylem compartment i [g h^-1]"],
        Wₓ(t)[1:nz], [description = "Water mass of root xylem compartment i [g]"],
        hₛ(t)[1:nz], [description = "Bulk soil matric potential in layer i [cm]"],
        dhₛ(t), [description = "Change in bulk soil matric potential in layer i [cm h^-1]"],
        Hₛ(t)[1:nz], [description = "Water potential of soil in layer i [cm]"],
        s(t)[1:nz], [description = "Water uptake from soil layer i [cm h^-1]"],
        r_rhiz(t)[1:nz], [description = "Root-soil interface radius in layer i [cm]"],
        rld(t)[1:nz], [description = "Root length density in layer i [cm cm^-3]"],
        ρ(t)[1:nz], [description = "Root-soil contact fraction in layer i [-]"],
        B(t)[1:nz], [description = "Root-soil interface conductance in layer i [cm h^-1]"],
        Hᵣₛ(t)[1:nz], [description = "Water potential at the root-soil interface in layer i [cm]"],
        Kᵣ(t)[1:nz], [description = "Radial root hydraulic conductivity in layer i [cm h^-1]"],
        uptake(t)[1:nz], [description = "Water uptake by the root from soil layer i [cm h^-1]"],
        H₀(t), [description = "Water potential at the root base [cm]"],
        LAI(t), [description = "Leaf area index [-]"],
        f_Ψ(t), [description = "Soil water stress factor [-]"],
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
        Tmin = Tmin, [description = "minimum temperature for growth [°C]"],
        Tmax = Tmax, [description = "maximum temperature for growth [°C]"],
        Topt = Topt, [description = "optimal temperature for growth [°C]"],
        v_max = v_max, [description = "maximum rate of vegetative development [h^-1]"],
        S_ref = S_ref, [description = "vegetative development stage at which reproductive development starts [-]"],
        k_s = k_s, [description = "sensitivity of development to vegetative development stage [-]"],
        k_Ψ = k_Ψ, [description = "sensitivity of development to water potential [-]"],
        Ψ_ref = Ψ_ref, [description = "reference water potential for LAI development [MPa]"],
        r_LAI = r_LAI, [description = "rate of change of LAI with respect to vegetative development stage [m^2 m^-2 h^-1]"],
        r_max = r_max, [description = "maximum rate of reproductive development [h^-1]"],
    )
    vars = @variables (
        f_T(t), [description = "Effect of temperature on development [-]"],
        f_R(t), [description = "Effect of vegetative development stage on reproductive development [-]"],
        f_Ψ(t), [description = "Effect of water potential on development [-]"],
        T(t), [description = "Air temperature [°C]"],
        Sᵥ(t), [description = "Vegetative development stage [-]"],
        dSᵥ(t), [description = "Rate of change of vegetative development stage [h^-1]"],
        LAI(t), [description = "Leaf area index [m^2 m^-2]"],
        dLAI(t), [description = "Rate of change of leaf area index [m^2 m^-2 h^-1]"],
        Ψ(t), [description = "Water potential [MPa]"],
        Sᵣ(t), [description = "Reproductive development stage [-]"],
        dSᵣ(t), [description = "Rate of change of reproductive development stage [h^-1]"],
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
        dz = dz, [description = "Layer width [cm]"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1 [cm h^-1]"],
        K_half(t), [description = "Hydraulic conductivity of connection [cm h^-1]"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1 [cm h^-1]"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2 [cm h^-1]"],
        Ψ_1(t), [description = "Total water potential of compartment 1 [cm]"],
        Ψ_2(t), [description = "Total water potential of compartment 2 [cm]"],
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