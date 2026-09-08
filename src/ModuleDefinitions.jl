@independent_variables t
D = Differential(t)

# Note on units:
# Everything is implicitly multiplied with 1 cm² soil surface area
# (cm^3 / cm^2) => cm

# # Modules
# ## Translated
function soil_module(; name, Ψ_m, α, n, Kₛ, l, θₛ, θᵣ, dz, z)
    ρ_w = 1.0 # Density of water [g cm^-3]
    g = 9.8 * 1.0e-5 # Gravitational acceleration [MPa cm^2 g^-1]
    Pₕ = ρ_w * g * z # Gravitational water potential [MPa]

    @parameters(
        α = α, [description = "van Genuchten shape parameter (related to inverse of air entry suction, > 0) [cm^-1]"],
        n = n, [description = "van Genuchten shape parameter (related to pore size distribution, > 1) [-]"],
        Kₛ = Kₛ, [description = "Saturated hydraulic conductivity [cm h^-1]"],
        l = l, [description = "Pore connectivity parameter [-]"],
        θₛ = θₛ, [description = "Saturated volumetric water content [-]"],
        θᵣ = θᵣ, [description = "Residual volumetric water content [-]"],
        dz = dz, [description = "Layer width [cm]"],
        Pₕ = Pₕ, [description = "Gravitational water potential [MPa]"],
    )
    @variables (
        Ψ(t), [description = "Total water potential [MPa]"], # alias `hT`
        Ψ_m(t) = Ψ_m, [description = "Matric water potential [MPa]"], # alias `h`
        C(t), [description = "Soil water capacitance [cm^-1]"],
        K(t), [description = "Hydraulic conductivity [cm h^-1]"], #eigenlijk moeten we dit zien als g per cm² per h, mits ρ_w = 1.0 g cm^-3
        θ(t), [description = "Volumetric water content [-]"],
        F(t), [description = "Water flux [cm h^-1]"], # alias `q` eigenlijk moeten we dit zien als g per cm² per h, mits ρ_w = 1.0 g cm^-3
        ΣF(t), [description = "Net water influx [cm h^-1]"], # alias `dq + s`  eigenlijk moeten we dit zien als g per cm² per h, mits ρ_w = 1.0 g cm^-3
        z(t) = z, [description = "Layer depth [cm]"],
    )
    eqs = [
        C ~ vanGenuchten_C(Ψ, θₛ, θᵣ, α, n),
        K ~ vanGenuchten_K(Ψ, θₛ, θᵣ, α, n, Kₛ, l),
        θ ~ vanGenuchten_θ(Ψ, θₛ, θᵣ, α, n),

        D(z) ~ 0,
        D(Ψ_m) ~ ( ΣF/dz ) / C,
        Ψ ~ Ψ_m + Pₕ # eq. 7 in paper
    ]

    system = ODESystem(eqs, t; name)
    return system
end

function rootuptake_module(; name, εₓ, rᵣ, dz, rld, Ψ)
    @constants ρ_w = 1.0 # density of water [g cm^-3]
    @parameters (
        εₓ = εₓ, [description = "Root xylem elastic modulus [MPa]"],
        rᵣ = rᵣ, [description = "Root radius [cm]"],
        dz = dz, [description = "Layer width [cm]"],
        rld = rld, [description = "Root length density [cm cm^-3]"], # the holy grail
    )
    @variables (        
        lᵣ(t), [description = "Root length in soil layer [cm]"],
        r_rhiz(t), [description = "Root-soil interface radius [cm]"],
        ρ(t), [description = "Root-soil contact fraction [-]"],
        B(t), [description = "Root-soil interface conductance [cm h^-1]"],
        Aᵣ(t), [description = "Normalized root surface area for root radial water transport [cm]"],
        Vᵣ(t), [description = "Normalized root cross-sectional area for root axial water transport [cm h^-1]"],
        Wₓ(t), [description = "Water mass of root xylem [g]"],
        dWₓ(t), [description = "Change in water mass of root xylem [g h^-1]"],
        Ψ(t) = Ψ, [description = "Water potential of root xylem [cm]"],
        dΨ(t), [description = "Change in water potential of root xylem [cm h^-1]"],
        ΣF(t), [description = "Net water influx [cm h^-1]"],
    )
    eqs = [
        lᵣ ~ rld * dz, # total root length
        r_rhiz ~ 1/sqrt(π * rld), # root-soil interface
        ρ ~ r_rhiz / rᵣ, # root-soil contact fraction
        B ~ 2*(ρ^2 - 1) / (1 - (0.53*ρ)^2 + 2*ρ^2*(log(0.53) + log(ρ))), # root-soil interface conductance
        Aᵣ ~ (2*π * rᵣ * lᵣ) / rᵣ, # normalized root surface area
        Vᵣ ~ 1 / dz * (lᵣ/dz * π * rᵣ^2), # normalized root cross-sectional area
        Wₓ ~ lᵣ * π * ρ_w * rᵣ^2, # representative "water mass" of the root xylem (g)
        dWₓ ~ ΣF, # change in water mass of the root xylem #! why not enforce D(Wₓ) ~ Wₓ (can be numerically unstable if params are wrong)
        dΨ ~ εₓ / Wₓ*dWₓ, # change in xylem water potential (cm h⁻¹)
        D(Ψ) ~ dΨ,
    ]
    
    system = ODESystem(eqs, t; name)
    return system
end

# ## TODO
#=
function rootbase()

    (
        Ψ₀(t), [description = "Water potential at the root base [cm]"],
        f_Ψ(t), [description = "Soil water stress factor [-]"],
        Tp(t), [description = "Transpiration rate [cm h^-1]"],
        Ψ_ref = Ψ_ref, [description = "Reference water potential [MPa]"],
        k_Ψ = k_Ψ, [description = "Water potential sensitivity coefficient [-]"],
        kc = kc, [description = "Crop coefficient for transpiration [-]"],
    )

    [
        # collar water potential (cm)
        Ψ₀ ~ max(Ψₓ[1] - Tp/(Kₓ[1] + tiny), hmin),
        #Ψ ~ Ψ₀ * 98.1e-6,
        f_Ψ ~ 0.5 *( 1 + tanh(k_Ψ * (Ψ₀ * 98.1e-6 - Ψ_ref))),
        Tp ~ ET0_inputfun(t) /10.0 * kc * (1 - exp(-0.45*(LAI))) * f_Ψ,
    ]

end

function phenology_module(; name, Tmin, Tmax, Topt, v_max, S_ref, k_s, k_Ψ, Ψ_ref, r_LAI, r_max)
    @parameters(
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
    @variables (
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
=#

# # Module connections
function soil_connection(; name, dz)
    @parameters(
        dz = dz, [description = "Layer width [cm]"],
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

function root_connection(; name, kₓ)
    @parameters(
        kₓ = kₓ, [description = "Intrinsic axial root hydraulic conductivity [h^-1]"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1 [cm h^-1]"],
        K_half(t), [description = "Hydraulic conductivity of connection [cm h^-1]"],
        Kₓ_1(t), [description = "Hydraulic xylem conductivity of compartment 1 [cm h^-1]"],
        Kₓ_2(t), [description = "Hydraulic xylem conductivity of compartment 2 [cm h^-1]"],
        Ψₓ_1(t), [description = "Water potential of root xylem compartment 1 [cm]"],
        Ψₓ_2(t), [description = "Water potential of root xylem compartment 2 [cm]"],
        Vᵣ_1(t), [description = "Normalized root volume of compartment 1 [cm h^-1]"],
        Vᵣ_2(t), [description = "Normalized root volume of compartment 2 [cm h^-1]"],
    )
    eqs = [
        F ~ K_half * (Ψₓ_2 - Ψₓ_1), 
        #! can simplifiy next 3 equations into 1 if we dont care about Kₓ
        K_half ~ Kₓ_1*Kₓ_2 * 2/(Kₓ_1 + Kₓ_2), 
        Kₓ_1 ~ kₓ * Vᵣ_1,
        Kₓ_2 ~ kₓ * Vᵣ_2,
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = [
        connection_MTK.Ψₓ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψₓ_2 ~ nb_node_MTK.Ψ,
        connection_MTK.Vᵣ_1 ~ node_MTK.Vᵣ,
        connection_MTK.Vᵣ_2 ~ nb_node_MTK.Vᵣ,
    ]
    return System(eqs, t; name), get_connection_eqset
end

function root_soil_connection(; name, kᵣ, hₛ, α, n, Kₛ, l)
    @parameters(
        kᵣ = kᵣ, [description = "Intrinsic radial root hydraulic conductivity [h^-1]"],
        hₛ = hₛ, [description = "Bulk soil matric potential [cm]"],
        α = α, [description = "van Genuchten shape parameter (related to inverse of air entry suction, > 0) [cm^-1]"],
        n = n, [description = "van Genuchten shape parameter (related to pore size distribution, > 1) [-]"],
        Kₛ = Kₛ, [description = "Saturated hydraulic conductivity [cm h^-1]"],
        l = l, [description = "Pore connectivity parameter [-]"],
    )
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1 [cm h^-1]"],
        Kᵣ(t), [description = "Radial root hydraulic conductivity (conductivity of connection) [cm h^-1]"],
        Ψₓ(t), [description = "Water potential of root xylem compartment [cm]"],
        Ψₛ(t), [description = "Water potential of soil [cm]"],
        Ψᵣₛ(t), [description = "Water potential at the root-soil interface [cm]"],
        Aᵣ(t), [description = "Normalized root surface area [cm]"],
        B(t), [description = "Root-soil interface conductance [cm h^-1]"],
        z(t), [description = "Layer depth [cm]"],
    )
    eqs = [
        F ~ (Ψᵣₛ - Ψₓ) * Kᵣ,
        Kᵣ ~ kᵣ * Aᵣ, # Hydraulic conductivities based on root dimensions and intrinsic k's (h⁻¹)
        kᵣ * Ψₓ + B * Ψₛ * ksoilfun(hₛ, Ψᵣₛ - z, α, n, Kₛ, l) ~ Ψᵣₛ * (B * ksoilfun(hₛ, Ψᵣₛ - z, α, n, Kₛ, l) + kᵣ),
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, original_order) = (
        original_order ?
        [
            connection_MTK.Ψₓ ~ node_MTK.Ψ,
            connection_MTK.Aᵣ ~ node_MTK.Aᵣ,
            connection_MTK.B ~ node_MTK.B,
            connection_MTK.Ψₛ ~ nb_node_MTK.Ψ,
            connection_MTK.z ~ nb_node_MTK.z,
        ] :
        [
            connection_MTK.Ψₓ ~ nb_node_MTK.Ψ,
            connection_MTK.Aᵣ ~ nb_node_MTK.Aᵣ,
            connection_MTK.B ~ nb_node_MTK.B,
            connection_MTK.Ψₛ ~ node_MTK.Ψ,
            connection_MTK.z ~ node_MTK.z,
        ]
    )
    
    return System(eqs, t; name), get_connection_eqset
end