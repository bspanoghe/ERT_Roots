@independent_variables t
D = Differential(t)
°C_to_K(T::Number) = T + 273.15 # temperature unit conversion

# Note on units:
# Everything is implicitly multiplied with 1 cm² soil surface area
# (cm^3 / cm^2) => cm

# # Modules
# ## Translated
function soil_module(; name, Ψ_m, α, n, Kₛ, l, θₛ, θᵣ, dz, z)
    ρ_w = 1.0 # Density of water [g cm^-3]
    g = MPa2cm(9.8e-5) # Gravitational acceleration [cm cm^2 g^-1]
    Pₕ = ρ_w * g * z # Gravitational water potential [cm]

    @parameters(
        α = α, [description = "van Genuchten shape parameter (related to inverse of air entry suction, > 0) [cm^-1]"],
        n = n, [description = "van Genuchten shape parameter (related to pore size distribution, > 1) [-]"],
        Kₛ = Kₛ, [description = "Saturated hydraulic conductivity [cm h^-1]"],
        l = l, [description = "Pore connectivity parameter [-]"],
        θₛ = θₛ, [description = "Saturated volumetric water content [-]"],
        θᵣ = θᵣ, [description = "Residual volumetric water content [-]"],
        dz = dz, [description = "Layer width [cm]"],
        Pₕ = Pₕ, [description = "Gravitational water potential [cm]"],
    )
    @variables (
        Ψ(t), [description = "Total water potential [cm]"], # alias `hT`
        Ψ_m(t) = Ψ_m, [description = "Matric water potential [cm]"], # alias `h`
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

    return ODESystem(eqs, t; name)
end

function rootuptake_module(; name, εₓ, rᵣ, dz, rld, Ψ)
    @constants ρ_w = 1.0 # density of water [g cm^-3]
    @parameters (
        εₓ = εₓ, [description = "Root xylem elastic modulus [cm]"],
        rᵣ = rᵣ, [description = "Root radius [cm]"],
        dz = dz, [description = "Layer width [cm]"],
        # rld = rld, [description = "Root length density [cm cm^-3]"],
    )
    @variables (
        rld(t) = rld, [description = "Root length density [cm cm^-3]"],
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
        B ~ 2*(ρ^2 - 1) / (1 - (0.53*ρ)^2 + 2*ρ^2*(log(0.53) + log(ρ))), # root-soil interface conductance #! 0.53?
        Aᵣ ~ (2*π * rᵣ * lᵣ) / rᵣ, # normalized root surface area
        Vᵣ ~ 1 / dz * (lᵣ/dz * π * rᵣ^2), # normalized root cross-sectional area
        Wₓ ~ lᵣ * π * ρ_w * rᵣ^2, # representative "water mass" of the root xylem (g)
        dWₓ ~ ΣF, # change in water mass of the root xylem
        dΨ ~ εₓ / Wₓ*dWₓ, # change in xylem water potential (cm h⁻¹)
        D(Wₓ) ~ dWₓ,
        D(Ψ) ~ dΨ,
    ]
    
    return ODESystem(eqs, t; name)
end

# goes together with rootuptake_module
function collar_module(; 
        name, k_Ψ_transp, kc, Tair, Tmin, Tmax, Topt, v_max, S_ref, k_s, k_Ψ_dev, Ψ_ref, r_LAI, r_max, Sᵥ, Sᵣ, LAI
    )
    @parameters(
        # root base
        Ψ_ref = Ψ_ref, [description = "Reference water potential [cm]"],
        k_Ψ_transp = k_Ψ_transp, [description = "Water potential sensitivity coefficient [-]"],
        kc = kc, [description = "Crop coefficient for transpiration [-]"],

        # phenology
        Tair = Tair, [description = "Air temperature [°C]"], #! changed to param
        Tmin = Tmin, [description = "minimum temperature for growth [°C]"],
        Tmax = Tmax, [description = "maximum temperature for growth [°C]"],
        Topt = Topt, [description = "optimal temperature for growth [°C]"],
        k_s = k_s, [description = "sensitivity of development to vegetative development stage [-]"],
        S_ref = S_ref, [description = "vegetative development stage at which reproductive development starts [-]"],
        k_Ψ_dev = k_Ψ_dev, [description = "sensitivity of development to water potential [-]"],
        v_max = v_max, [description = "maximum rate of vegetative development [h^-1]"],
        r_max = r_max, [description = "maximum rate of reproductive development [h^-1]"],
        r_LAI = r_LAI, [description = "rate of change of LAI with respect to vegetative development stage [m^2 m^-2 h^-1]"],
    )
    @variables (
        # root base
        Ψ(t), [description = "Water potential of root xylem [cm]"], # defined in rootuptake_module
        f_Ψ_transp(t), [description = "Soil water stress factor [-]"],
        f_LAI(t), [description = "Effect of LAI on transpiration [-]"], #! correct?
        Tp(t), [description = "Transpiration rate [cm h^-1]"],

        # phenology
        f_T(t), [description = "Effect of temperature on development [-]"],
        f_R(t), [description = "Effect of vegetative development stage on reproductive development [-]"],
        f_Ψ_dev(t), [description = "Effect of water potential on development [-]"],
        dSᵥ(t), [description = "Rate of change of vegetative development stage [h^-1]"],
        dSᵣ(t), [description = "Rate of change of reproductive development stage [h^-1]"],
        dLAI(t), [description = "Rate of change of leaf area index [m^2 m^-2 h^-1]"],
        Sᵥ(t) = Sᵥ, [description = "Vegetative development stage [-]"],
        Sᵣ(t) = Sᵣ, [description = "Reproductive development stage [-]"],
        LAI(t) = LAI, [description = "Leaf area index [m^2 m^-2]"],
    )
    eqs = [
        f_Ψ_transp ~ 0.5 *( 1 + tanh(k_Ψ_transp * (Ψ - Ψ_ref))),
        f_LAI ~ 1 - exp(-0.45*(LAI)),
        Tp ~ ET0_inputfun(t) * kc * f_LAI * f_Ψ_transp,

        f_T ~ (Tmax - Tair)/(Tmax - Topt) * ((Tair - Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)),
        f_R ~ smoothstep(k_s * (S_ref - Sᵥ)),
        f_Ψ_dev ~ smoothstep(k_Ψ_dev * (Ψ - Ψ_ref)),
        dSᵥ ~ v_max * f_T * f_R * smoothstep(10 * (t - 1000)), #!
        dSᵣ ~ (1-f_R) * r_max * f_T ,
        dLAI ~ dSᵥ * r_LAI * Sᵥ * (S_ref - Sᵥ)/S_ref * f_Ψ_dev,
        D(Sᵥ) ~ dSᵥ,
        D(Sᵣ) ~ dSᵣ, 
        D(LAI) ~ dLAI,
    ]
    system = ODESystem(eqs, t; name)

    return system
end

function Ψ_air_module_cm(; name, T)
    @variables (
        Ψ(t), [description = "Total water potential [cm]"],
        W_r(t), [description = "Relative water content [g/g]"],
    )
    @parameters T = T [description = "Temperature [°C]"]
    @constants (
        R = MPa2cm(8.314), [description = "Ideal gas constant [cm cm^3 K^-1 mol^-1]"],
        V_w = 18, [description = "Molar volume of water [cm^3/mol]"],
    )

    eqs = [Ψ ~ R * °C_to_K(T) / V_w * log(W_r)] # Spanner equation (see e.g. https://academic.oup.com/insilicoplants/article/4/1/diab038/6510844)

    return System(eqs, t; name)
end

function Ψ_soil_module_cm(; name)
    @variables (
        Ψ(t), [description = "Total water potential [cm]"],
        W_r(t), [description = "Relative water content [g g^-1]"],
    )

    eqs = [Ψ ~ MPa2cm(soilfunc(W_r))]

    return System(eqs, t; name)
end
soilfunc(W_r; a = 3.5, b = 5.5) = -(a / W_r) * exp(-b * W_r) # empirical equation for soil water potential
# default values fit to loam soil data from Chen et al. (1997)
# link: https://doi.org/10.1093/treephys/17.12.797

# # Module connections
function soil_connection(; name, dz)
    @parameters(
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

function collar_air_connection(; name, original_order)
    @variables (
        F(t), [description = "Water flux from compartment 2 (air) to compartment 1 (collar) [cm h^-1]"],
        Tp(t), [description = "Transpiration rate [cm h^-1]"],
    )

    polarity = original_order ? 1 : -1

    eqs = [
        F ~ polarity * Tp,
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK, original_order) = (
        original_order ?    
        [connection_MTK.Tp ~ node_MTK.Tp] :
        [connection_MTK.Tp ~ nb_node_MTK.Tp]
    ) 
    return System(eqs, t; name), get_connection_eqset
end