using Pkg; Pkg.activate("./scripts")
using ModelingToolkit, OrdinaryDiffEq, PlantModules, PlantModules.PlantGraphs, Plots

# # Structure

struct Air <: Node end
struct Soil{T} <: Node
    z::T
end
struct Drainage <: Node end

struct Root{T} <: Node
    z::T
    rld::T
end

const dz = 1.0
n_layers = 3
add_roots = true
soil_graph = sum([Soil(-i*dz) for i in 1:n_layers])

if !add_roots
    graphs = [Air(), soil_graph, Drainage()]

    intergraph_connections = [
        (1, 2) => (:Air, getnodes(soil_graph)[1]),
        (2, 3) => (getnodes(soil_graph)[end], :Drainage),
    ]
else
    root_graph = sum([Root(-i*dz, 1.0) for i in 1:n_layers])
    graphs = [Air(), soil_graph, Drainage(), root_graph]

    if n_layers == 1
        is_root_soil_connected(root, soil) = true
    else
        is_root_soil_connected(root, soil) = data(root).z == data(soil).z
    end

    intergraph_connections = [
        (1, 2) => (:Air, getnodes(soil_graph)[1]),
        (2, 3) => (getnodes(soil_graph)[end], :Drainage),
        (2, 4) => is_root_soil_connected,
    ]
end

plantstructure = PlantStructure(graphs, intergraph_connections)
plotstructure(plantstructure)


# # Function

# ## Functional modules

# ### Putting the fun in functions

# Calculate volumentric water content θ from pressure head h and soil hydraulic parameters
function vanGenuchten_θ(h, θ_s, θ_r, α, n)
    m = 1 - 1 / n 
    h_eff = -h # -min(h, -1e-6)
    # Volumetric moisture content (θ)
    θ = (θ_s - θ_r) / (1 + (α * h_eff)^n)^m + θ_r
    return θ
end

# Calculate hydraulic conductivity K from pressure head h and soil hydraulic parameters
function vanGenuchten_K(h, θ_s, θ_r, α, n, K_s, l)
    m = 1 - 1 / n 
    h_eff = -h # -min(h, -1e-6)
    # Volumetric moisture content (θ)
    θ = (θ_s - θ_r) / (1 + (α * h_eff)^n)^m + θ_r
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
    h_eff = -h # -min(h, -1e-6)
    C = (θ_s - θ_r)* n * m * α * (α * h_eff)^(n - 1) * 
         (1 + (α * h_eff)^n)^(1/n - 2)
    return C
end

# @register_symbolic vanGenuchten_θ(h, θ_s, θ_r, α, n)
# @register_symbolic vanGenuchten_K(h, θ_s, θ_r, α, n, K_s, l)
# @register_symbolic vanGenuchten_C(h, θ_s, θ_r, α, n)

function ksoilfun(hsoil, hint, α, n, Ks, l)
    (fluxmpfunction(hsoil, α, n, Ks, l) - fluxmpfunction(hint, α, n, Ks, l)) / (hsoil - hint)# + eps())
end
fluxmpfunction(hub, α, n, Ks, l; hlb = -1000.0) = integrate(h -> Khfunc(h, α, n, Ks, l), hlb, hub)
function Khfunc(h, α, n, Ks, l)
    m = 1 - 1 / n
    α_h = max(zero(h), -α*h) #! curse you domain errors
    q = Ks * (1 - α_h^(n - 1) * (1 + α_h^n)^(-m))^2 *
        (1 + α_h^n)^(-l * m)
    return q
end
@register_symbolic ksoilfun(hsoil, hint, α, n, Ks, l)

# placeholder to minimize amount of dependencies
function integrate(f, lb, ub; n = 1000)
    int_range = range(lb, ub, length = n)
    Δx = step(int_range)
    sum(f.(int_range)) * Δx
end




# quick tests
finesoil = (θₛ = 0.43, θᵣ = 0.078, α = 0.0083, n = 1.2539, Kₛ = 2.272 / (24), l = 0.5)
plot(h -> vanGenuchten_θ(h, values(finesoil)[1:4]...), xlims = (-1000.0, 100.0))
plot(h -> vanGenuchten_K(h, values(finesoil)...), xlims = (-1000.0, 100.0))
plot(h -> vanGenuchten_C(h, values(finesoil)[1:4]...), xlims = (-1000.0, 100.0))

plot(h -> Khfunc(h, values(finesoil)[3:end]...), xlims = (-1000, 100))
plot(h -> ksoilfun(h, -10.0, values(finesoil)[3:end]...), xlims = (-1000, 100), ylims = (0.0, 0.05))

# ### Modules

include("../src/ModuleDefinitions.jl")

# ## Coupling

module_coupling = Dict(
    :Soil => [soil_module],
    :Air => [environmental_module, Ψ_air_module],
    :Drainage => [environmental_module, Ψ_soil_module],
    :Root => [rootuptake_module],
);
connecting_modules = Dict(
    (:Air, :Soil) => constant_hydraulic_connection,
    (:Soil, :Soil) => soil_connection,
    (:Soil, :Drainage) => constant_hydraulic_connection,
    (:Root, :Root) => root_connection,
    (:Root, :Soil) => root_soil_connection
);

plantcoupling = PlantCoupling(; module_coupling, connecting_modules);

# ## Parameters
cm2MPa(x) = 98.1e-6 * x
MPa2cm(x) = 1/98.1e-6 * x

default_changes = Dict(
    Pair.(keys(finesoil), values(finesoil))..., 
    :dz => dz,
    :W_max => 1e5,
    :z => NaN, # assigned in nodes
    :Ψ_m => MPa2cm(-3.0),
    :Ψ => MPa2cm(-1.0),
    :εₓ => MPa2cm(10.0),
    :rᵣ => 0.05,
    :dz => dz, 
    :rld => NaN, # assigned in nodes
    :kₓ => 0.1,
    :kᵣ => 0.1,
    :hₛ => MPa2cm(-3.0)
);
module_defaults = Dict(
    :Air => Dict(:W_r => 0.9),
    :Soil => Dict(:Kₛ => 100.0),#, :Kₛ => 10.0, :l => 2.0),
    :Drainage => Dict(:W_r => 0.1)
);
connection_values = Dict(
    (:Air, :Soil) => Dict(:K => 1e-5),
    (:Soil, :Drainage) => Dict(:K => 1e-5)
);
plantparams = PlantParameters(; default_changes, module_defaults, connection_values);

system = generate_system(plantstructure, plantcoupling, plantparams)
time_span = (0.0, 48.0);
prob = ODEProblem(system, [], time_span, sparse = true);
sol = solve(prob, FBDF());
plotgraph(sol, plantstructure, varname = :θ, structmod = :Soil)
plotgraph(sol, plantstructure, varname = :Wₓ, structmod = :Root, ylims = (0, 0.1))
plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Soil, :Drainage, :Air, :Root], ylims = (-100, 0))
plotgraph(sol, plantstructure, varname = :Ψ, structmod = [:Soil, :Drainage, :Air, :Root])

plotgraph(sol, plantstructure, varname = :ΣF, structmod = [:Air, :Drainage])
plotgraph(sol, plantstructure, varname = :W, structmod = :Drainage)