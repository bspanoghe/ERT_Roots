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
include("../src/ModuleDefinitions.jl")

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