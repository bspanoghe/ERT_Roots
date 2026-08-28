using Pkg; Pkg.activate("./scripts")
using ModelingToolkit, OrdinaryDiffEq, PlantModules, PlantModules.PlantGraphs, Plots

# # Structure

struct Air <: Node end
struct Soil{T} <: Node
    z::T
end
struct Drainage <: Node end

const dz = 1.0
soil_graph = sum([Soil(-i*dz) for i in 1:3])
full_graph = Air() + soil_graph + Drainage()
plantstructure = PlantStructure(full_graph)
plotstructure(plantstructure)

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
finesoil = (θₛ = 0.43, θᵣ = 0.078, α = 0.0083, n = 1.2539, Kₛ = 2.272 / (24), l = 0.5)
plot(h -> vanGenuchten_θ(h, values(finesoil)[1:4]...), xlims = (-1000.0, 0.0))
plot(h -> vanGenuchten_K(h, values(finesoil)...), xlims = (-1000.0, 0.0))
plot(h -> vanGenuchten_C(h, values(finesoil)[1:4]...), xlims = (-1000.0, 0.0))

# ### Modules
@independent_variables t
D = Differential(t)

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

# ### Module connections
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

# ## Coupling

module_coupling = Dict(
    :Soil => [soil_module],
    :Air => [environmental_module, Ψ_air_module],
    :Drainage => [environmental_module, Ψ_soil_module]
);
connecting_modules = Dict(
    (:Air, :Soil) => constant_hydraulic_connection,
    (:Soil, :Soil) => soil_connection,
    (:Soil, :Drainage) => constant_hydraulic_connection
)

plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ## Parameters

default_changes = Dict(
    Pair.(keys(finesoil), values(finesoil))..., 
    :dz => dz,
    :W_max => 1e5,
    :z => NaN, # must be assigned in nodes!
    :Ψ_m => -1.0
)
module_defaults = Dict(
    :Air => Dict(:W_r => 0.9),
    :Drainage => Dict(:W_r => 0.01)
)
connection_values = Dict(
    (:Air, :Soil) => Dict(:K => 1e-5),
    (:Soil, :Drainage) => Dict(:K => 1e-3)
)
plantparams = PlantParameters(; default_changes, module_defaults, connection_values);

system = generate_system(plantstructure, plantcoupling, plantparams)
time_span = (0.0, 48.0);
prob = ODEProblem(system, [], time_span, sparse = true);
sol = solve(prob, FBDF());
plotgraph(sol, plantstructure, varname = :θ, structmod = :Soil)
plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Soil, :Drainage])
plotgraph(sol, plantstructure, varname = :ΣF, structmod = [:Air, :Drainage])
plotgraph(sol, plantstructure, varname = :W, structmod = :Drainage)