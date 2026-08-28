using Pkg; Pkg.activate("./scripts")
using ModelingToolkit, OrdinaryDiffEq, PlantModules, PlantModules.PlantGraphs, Plots

# # Structure

struct Air <: Node end
struct Soil{T} <: Node
    z::T
end
struct Drainage <: Node end
soil_graph = Air() + Soil(-1.0) + Soil(-2.0) + Soil(-3.0) + Drainage()
plantstructure = PlantStructure(soil_graph)

# # Function

# ## Functional modules

# ### Putting the fun in functions

# Calculate volumentric water content θ from pressure head h and soil hydraulic parameters
function vanGenuchten_θ(h, θ_s, θ_r, α, n)
    m = 1 - 1 / n 
    h_eff = min(h, -1e-6)
    # Volumetric moisture content (θ)
    θ = (θ_s - θ_r) / (1 + (α * abs(h_eff))^n)^m + θ_r
    return θ
end

# Calculate hydraulic conductivity K from pressure head h and soil hydraulic parameters
function vanGenuchten_K(h, θ_s, θ_r, α, n, K_s, l)
    m = 1 - 1 / n 
    h_eff = min(h, -1e-6)
    # Volumetric moisture content (θ)
    θ = (θ_s - θ_r) / (1 + (α * abs(h_eff))^n)^m + θ_r
    # Effective saturation (Se)
    Se = (θ - θ_r) / (θ_s - θ_r)
    # Hydraulic conductivity (K)
    K = K_s * Se^l * (1 - (1 - Se^(1 / m))^m)^2
    return  K
end

# Calculate specific soil water capacitance C from pressure head h and soil hydraulic parameters
function vanGenuchten_C(h, θ_s, θ_r, α, n)
    # Specific moisture storage (C)
    m = 1 - 1 / n 
    h_eff = min(h, -1e-6)
    C = (θ_s - θ_r)* n * m * α * (α * abs(h_eff))^(n - 1) * 
         (1 + (α * abs(h_eff))^n)^(1/n - 2)
    return C
end

@register_symbolic vanGenuchten_θ(h, θ_s, θ_r, α, n)
@register_symbolic vanGenuchten_K(h, θ_s, θ_r, α, n, K_s, l)
@register_symbolic vanGenuchten_C(h, θ_s, θ_r, α, n)

# quick tests

θ_s = 0.43
θ_r = 0.078
α = 0.0083
n = 1.2539
K_s = 2.272 / (24)
l = 0.5

plot(h -> vanGenuchten_θ(h, θ_s, θ_r, α, n), xlims = (-3000.0, 0.0))
plot(h -> vanGenuchten_K(h, θ_s, θ_r, α, n, K_s, l), xlims = (-3000.0, 0.0))
plot(h -> vanGenuchten_C(h, θ_s, θ_r, α, n), xlims = (-3000.0, 0.0))

# ### Actual modules

@independent_variables t
D = Differential(t)

function soil_module(; name, α, n, Kₛ, l, θₛ, θᵣ, dz)
    params = @parameters(
        α = α, [description = ""],
        n = n, [description = ""],
        Kₛ = Kₛ, [description = ""],
        l = l, [description = ""],
        θₛ = θₛ, [description = ""],
        θᵣ = θᵣ, [description = ""],
        dz = dz, [description = ""],
    )
    vars = @variables (
        Ψ(t), [description = "Total water potential [?]"], # alias `hT`
        Ψ_m(t), [description = "Matrix water potential [?]"], # alias `h`
        C(t), [description = "Soil water capacitance [?]"],
        K(t), [description = "Hydraulic conductivity [?]"],
        θ(t), [description = "Volumetric water content [?]"],
        s(t), [description = "Root water uptake sink term (?) [?]"],
        F(t), [description = "Water flux"], # alias `q`
        dF(t), [description = "Change in water flux"], # alias `dq`
    )
    eqs = [
        C ~ vanGenuchten_C(Ψ, θₛ, θᵣ, α, n),
        K ~ vanGenuchten_K(Ψ, θₛ, θᵣ, α, n, Kₛ, l),
        θ ~ vanGenuchten_θ(Ψ, θₛ, θᵣ, α, n),

        # time derivative: C * dh/dt = -dq - s  (signs depend on conventions)
        D(Ψ[i]) ~ ( dF/ dz - s/dz ) / C,
        [hT[i] ~ z[i] + h[i] for i in 1:nz]...,        # eq. 7 in paper
    ]

    system = ODESystem(eqs, t; name)
    return system
end

function soil_connection(; name)
    @variables (
        F(t), [description = "Water flux from compartment 2 to compartment 1"], #, unit = u"g / hr"],
        dF(t), [description = "Change in water flux from compartment 2 to compartment 1"], #, unit = u"g / hr"],
        K_half(t), [description = "Hydraulic conductivity of connection"], #, unit = u"g / hr / MPa"],
        K_1(t), [description = "Hydraulic conductivity of compartment 1"], #, unit = u"g / hr / MPa"],
        K_2(t), [description = "Hydraulic conductivity of compartment 2"], #, unit = u"g / hr / MPa"],
        Ψ_1(t), [description = "Total water potential of compartment 1"], #, unit = u"MPa"],
        Ψ_2(t), [description = "Total water potential of compartment 2"], #, unit = u"MPa"],
    )

    eqs = [
        F ~ K_half * ( (Ψ_2-Ψ_1) / dz -  dz/dz ),
        K_half ~ 2 / (1/K_1 + 1/K_2),
        D(F) ~ dF
    ]

    get_connection_eqset(node_MTK, nb_node_MTK, connection_MTK) = [
        connection_MTK.Ψ_1 ~ node_MTK.Ψ,
        connection_MTK.Ψ_2 ~ nb_node_MTK.Ψ,
        connection_MTK.K_1 ~ node_MTK.K,
        connection_MTK.K_2 ~ nb_node_MTK.K,
        connection_MTK.dF ~ node_MTK.dF, #! double check sign
        connection_MTK.dF ~ -nb_node_MTK.dF, #! double check sign
    ]
    return System(eqs, t; name), get_connection_eqset
end

# ## Parameters
