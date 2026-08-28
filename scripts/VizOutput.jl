# Plot the soil water potential across soil layers
plot(sol_ode, idxs=soil_fine.h, xlabel="Time", ylabel="h(t)", legend=:none)
savefig("./img/soilh.png")

# Plot the collar water potential 
plot(sol_ode, idxs=wup.H₀, xlabel="Time", ylabel="Collar Ψ (cm)", legend=:none)
savefig("./img/Hcollar.png")

# Plot the total water uptake and Potential Transpiration
summed_s = [sum(x) for x in sol_ode[wup.s]]
plot(0:1:3600, Transpiration.(0:1:3600), ylabel="(Potential) Water uptake (cm h⁻¹)")
plot!(0:1:3600, summed_s, xlabel="Time (h)", legend=:none)
savefig("./img/wateruptake.png")

z = .-z
t = sol_ode.t
H = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    # Fill row i (depth i) across time
    # If you use saved times:
    H[i, :] = [sol_ode(tt, idxs=soil_sl.h[i]) for tt in t]
    # Alternatively, if t == sol_ode.t, you can use sol_ode[soil.h[i]] depending on your setup;
    # the idxs route is the most general and works with interpolation.
end


heatmap(
    t, z, H;
    c = :viridis,              # viridis colormap
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Soil water potential (h)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)

thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_ode(tt, idxs=soil_sl.θ[i]) for tt in t]
end

heatmap(
    t, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison
    colorbar = false, 
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Volumetric water content (θ)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)


## Scenario 2: d_e = 30 io 65
z = .-z
t = sol_ode_2.t
thetas_2 = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas_2[i, :] = [sol_ode_2(tt, idxs=soil_fine.θ[i]) for tt in t]
end

heatmap(
    t, z, thetas_2;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Volumetric water content (θ)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)

## Scenario 3: t_e = 1800, t_m = 1200
z = .-z
t = sol_ode_3.t
thetas_3 = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas_3[i, :] = [sol_ode_3(tt, idxs=soil_fine.θ[i]) for tt in t]
end

heatmap(
    t, z, thetas_3;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Volumetric water content (θ)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)


## Scenario 4: t_e = 3600, t_m = 2400
z = .-z
t = sol_ode_4.t
thetas_4 = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas_4[i, :] = [sol_ode_4(tt, idxs=soil_fine.θ[i]) for tt in t]
end

heatmap(
    t, z, thetas_4;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Volumetric water content (θ)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)



##RLD
rld_1 = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld_1[i, :] = [sol_ode(tt, idxs=wup.rld[i]) for tt in t]
end

heatmap(
    t, z, rld_1;
    c = cgrad(:viridis, rev=true),              # viridis colormap
    clims = (0.0, 1.3), # consistent color scale for comparison
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Root length density (cm cm⁻³)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)


rld_3 = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld_3[i, :] = [sol_ode_3(tt, idxs=wup.rld[i]) for tt in t]
end

heatmap(
    t, z, rld_3;
    c = cgrad(:viridis, rev=true),              # viridis colormap
    clims = (0.0, 1.3), # consistent color scale for comparison
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Root length density (cm cm⁻³)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)

rld_4 = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld_4[i, :] = [sol_ode_4(tt, idxs=wup.rld[i]) for tt in t]
end

heatmap(
    t, z, rld_4;
    c = cgrad(:viridis, rev=true),              # viridis colormap
    clims = (0.0, 1.3), # consistent color scale for comparison
    xlabel = "Time",
    ylabel = "Depth (m)",
    colorbar_title = "Root length density (cm cm⁻³)",
    yflip = true,              # depth downward
    right_margin = 5Plots.mm,  # nicer colorbar spacing
)



##grouped figure

L1 = heatmap(
    title = "Volumetric water content (-)",
    t, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

L2 = heatmap(
    t, z, thetas_3;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison 
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true, 
    tickfontsize = 10, 
        legend_font_pointsize = 10,
    guidefontsize = 12              # depth downward
)
L3 = heatmap(
    t, z, thetas_4;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.1, 0.43), # consistent color scale for comparison
    xlabel = "Time (h)",
    ylabel = "Soil depth (cm)",
    yflip = true,
    tickfontsize = 10, 
        legend_font_pointsize = 10,
    guidefontsize = 12               # depth downward
)

R1 = heatmap(
    title = "Root length density (cm cm⁻³)",
    t, z, rld_1;
    c = cgrad(:viridis, rev=true),              # viridis colormap
    clims = (0.0, 1.3), # consistent color scale for comparison
    xlabel = "",
    ylabel = "",
    yflip = true, 
    tickfontsize = 10,
        legend_font_pointsize = 10,
    guidefontsize = 12               # depth downward
)

R2 = heatmap(
    t, z, rld_3;
    c = cgrad(:viridis, rev=true),              # viridis colormap
    clims = (0.0, 1.3), # consistent color scale for comparison
    xlabel = "",
    ylabel = "",
    yflip = true,  
    tickfontsize = 10,    
        legend_font_pointsize = 10,
    guidefontsize = 12          # depth downward# nicer colorbar spacing
)
R3 = heatmap(
    t, z, rld_4;
    c = cgrad(:viridis, rev=true),              # viridis colormap
    clims = (0.0, 1.3), # consistent color scale for comparison
    xlabel = "Time (h)",
    ylabel = nothing,
    yflip = true,  
    tickfontsize = 10, 
    legend_font_pointsize = 10,
    guidefontsize = 12          # depth downward
)


l = @layout [a{0.5w} b{0.5w}; c{0.5w} d{0.5w}; e{0.5w} f{0.5w}]

p =plot(L1, R1, L2, R2, L3, R3; layout = l, 
 tickfontsize = 10, size = (900, 600),
 guidefontsize = 12, legend_font_pointsize = 10, 
 legend_foreground_color = nothing,
left_margin = (4.0, :mm))

savefig(p, "./img/theta_rld.png")