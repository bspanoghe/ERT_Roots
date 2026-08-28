
#z = .-z
t = sol_ref.t
thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_ref(tt, idxs=soil_sl.θ[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref.t)
L1 = heatmap(
    Timestamp, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 0.3), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

t = sol_sh.t
thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_sh(tt, idxs=soil_sl.θ[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_sh.t)
R1 = heatmap(
    Timestamp, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 0.3), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)


#z = .-z
t = sol_ref_root.t
thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_ref_root(tt, idxs=soil_sl.θ[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref_root.t)
L2 = heatmap(
    Timestamp, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 0.3), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

t = sol_sh_root.t
thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_sh_root(tt, idxs=soil_sl.θ[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_sh_root.t)
R2 = heatmap(

    Timestamp, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 0.3), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)


#z = .-z
t = sol_ref_sto.t
thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_ref_sto(tt, idxs=soil_sl.θ[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref_sto.t)
L3 = heatmap(
    Timestamp, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 0.3), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

t = sol_sh_sto.t
thetas = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    thetas[i, :] = [sol_sh_sto(tt, idxs=soil_sl.θ[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_sh_sto.t)
R3 = heatmap(
    Timestamp, z, thetas;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 0.3), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

l = @layout [a{0.5w} b{0.5w}; c{0.5w} d{0.5w}; e{0.5w} f{0.5w}]

tm_ticks = round.(Timestamp, Month(1)) |> unique
p =plot(L1, R1, L2, R2, L3, R3; layout = l, 
 tickfontsize = 10, size = (900, 600),
 guidefontsize = 12, legend_font_pointsize = 10, 
 legend_foreground_color = nothing, 
 xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")),
left_margin = (4.0, :mm))

savefig(p, "./img/_Divagro_theta_rld.png")


df_Ψ = DataFrame(
        ref_Ψ = sol_ref[phenology.Ψ],
        ref_root_Ψ = sol_ref_root[phenology.Ψ],
        ref_sto_Ψ = sol_ref_sto[phenology.Ψ],
        sh_Ψ = sol_sh[phenology.Ψ],
        sh_root_Ψ = sol_sh_root[phenology.Ψ],
        sh_sto_Ψ = sol_sh_sto[phenology.Ψ],
        Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref_sto.t))

L1 = plot(df_Ψ.Timestamp, df_Ψ.ref_Ψ,
    xlabel = "",
    ylabel = "Water potential (MPa)",   
    tickfontsize = 10, color = :black,
    ylim = (-1.5, 0.0), label = :none
    ) 
L2 = plot(df_Ψ.Timestamp, df_Ψ.ref_root_Ψ,
    xlabel = "",
    ylabel = "Water potential (MPa)",   
    tickfontsize = 10, color = :black,
    ylim = (-1.5, 0.0), label = :none
    ) 

L3 = plot(df_Ψ.Timestamp, df_Ψ.ref_sto_Ψ,
    xlabel = "Time (h)",
    ylabel = "Water potential (MPa)",   
    tickfontsize = 10, color = :black,
    ylim = (-1.5, 0.0), label = :none
    )

R1 = plot(df_Ψ.Timestamp, df_Ψ.sh_Ψ,
    xlabel = "",
    ylabel = "Water potential (MPa)",   
    tickfontsize = 10, color = :black,
    ylim = (-1.5, 0.0), label = :none
    ) 
R2 = plot(df_Ψ.Timestamp, df_Ψ.sh_root_Ψ,
    xlabel = "",
    ylabel = "Water potential (MPa)",   
    tickfontsize = 10, color = :black,
    ylim = (-1.5, 0.0), label = :none
    ) 

R3 = plot(df_Ψ.Timestamp, df_Ψ.sh_sto_Ψ,
    xlabel = "Time (h)",
    ylabel = "Water potential (MPa)",   
    tickfontsize = 10, color = :black,
    ylim = (-1.5, 0.0), label = :none
    ) 

 p =plot(L1, R1, L2, R2, L3, R3; layout = l, 
 tickfontsize = 10, size = (900, 600),
 guidefontsize = 12, legend_font_pointsize = 10, 
 legend_foreground_color = nothing, 
 xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")),
left_margin = (4.0, :mm))

savefig(p, "./img/_Divagro_waterpotential.png")





t = sol_ref.t
rld = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld[i, :] = [sol_ref(tt, idxs=wup.rld[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref.t)
L1 = heatmap(
    Timestamp, z, rld;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 1.5), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

t = sol_sh.t
rld = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld[i, :] = [sol_sh(tt, idxs=wup.rld[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_sh.t)
R1 = heatmap(
    Timestamp, z, rld;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 1.5), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)


#z = .-z
t = sol_ref_root.t
rld = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld[i, :] = [sol_ref_root(tt, idxs=wup.rld[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref_root.t)
L2 = heatmap(
    Timestamp, z, rld;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 1.5), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

t = sol_sh_root.t
rld = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld[i, :] = [sol_sh_root(tt, idxs=wup.rld[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_sh_root.t)
R2 = heatmap(

    Timestamp, z, rld;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 1.5), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)


#z = .-z
t = sol_ref_sto.t
rld = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld[i, :] = [sol_ref_sto(tt, idxs=wup.rld[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref_sto.t)
L3 = heatmap(
    Timestamp, z, rld;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 1.5), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

t = sol_sh_sto.t
rld = Array{Float64}(undef, length(z), length(t))
for (i, zi) in enumerate(z)
    rld[i, :] = [sol_sh_sto(tt, idxs=wup.rld[i]) for tt in t]
end
Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_sh_sto.t)
R3 = heatmap(
    Timestamp, z, rld;
    c = cgrad(:thermal, rev=true),              # viridis colormap
    clims = (0.0, 1.5), # consistent color scale for comparison           # <-- puts legend on top (works in GR; if not, see note below)
    xlabel = "",
    ylabel = "Soil depth (cm)",
    yflip = true,   
    tickfontsize = 10,
    legend_font_pointsize = 10,
    guidefontsize = 12             # depth downward
)

l = @layout [a{0.5w} b{0.5w}; c{0.5w} d{0.5w}; e{0.5w} f{0.5w}]

tm_ticks = round.(Timestamp, Month(1)) |> unique
p =plot(L1, R1, L2, R2, L3, R3; layout = l, 
 tickfontsize = 10, size = (900, 600),
 guidefontsize = 12, legend_font_pointsize = 10, 
 legend_foreground_color = nothing, 
 xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")),
left_margin = (4.0, :mm))

savefig(p, "./img/_Divagro_rld.png")


Ψ_fun(Ψ, Ψ_ref, k_Ψ) = 0.5 *( 1 + tanh(k_Ψ * (Ψ - Ψ_ref)))

plot(collect(-1.5:0.01:0.0), Ψ_fun.(collect(-1.5:0.01:0.0), -0.75, 4.0), linewidth = 2, label = "reference",
xlabel = "Water potential (MPa)", ylabel = "f(Ψ)",
tickfontsize = 12,
 guidefontsize = 14, legend_font_pointsize = 12, 
 legend_foreground_color = nothing)
plot!(collect(-1.5:0.01:0.0), Ψ_fun.(collect(-1.5:0.01:0.0), -0.5, 5.0), linewidth = 2, label = "water-saving")

savefig("./img/_Divagro_sto.png")




df_Ψ = DataFrame(
        ref_Tp = sol_ref[wup.Tp],
        ref_root_Tp = sol_ref_root[wup.Tp],
        ref_sto_Tp = sol_ref_sto[wup.Tp],
        sh_Tp = sol_sh[wup.Tp],
        sh_root_Tp = sol_sh_root[wup.Tp],
        sh_sto_Tp = sol_sh_sto[wup.Tp],
        Timestamp = Ref_df.Timestamp[1] .+ Hour.(sol_ref_sto.t))

L1 = plot(df_Ψ.Timestamp, df_Ψ.ref_Tp,
    xlabel = "",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    ) 
L2 = plot(df_Ψ.Timestamp, df_Ψ.ref_root_Tp,
    xlabel = "",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    ) 

L3 = plot(df_Ψ.Timestamp, df_Ψ.ref_sto_Tp,
    xlabel = "Time (h)",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    )

R1 = plot(df_Ψ.Timestamp, df_Ψ.sh_Tp,
    xlabel = "",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    )

R1 = plot(df_Ψ.Timestamp, df_Ψ.sh_Tp,
    xlabel = "",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    ) 
R2 = plot(df_Ψ.Timestamp, df_Ψ.sh_root_Tp,
    xlabel = "",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    ) 

R3 = plot(df_Ψ.Timestamp, df_Ψ.sh_sto_Tp,
    xlabel = "Time (h)",
    ylabel = "Transp (cm h⁻¹)",   
    tickfontsize = 10, color = :black,
    ylim = (0, 0.1), label = :none
    ) 

 p =plot(L1, R1, L2, R2, L3, R3; layout = l, 
 tickfontsize = 10, size = (900, 600),
 guidefontsize = 12, legend_font_pointsize = 10, 
 legend_foreground_color = nothing, 
 xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")),
left_margin = (4.0, :mm))

savefig(p, "./img/_Divagro_transpiration.png")


t = collect(800:1:4200)
Timestamp = Ref_df.Timestamp[1] .+ Hour.(t)
tm_ticks = round.(Timestamp, Month(1)) |> unique

p1 = plot(Timestamp, Precip_ref.(t), linewidth = 2, color = :blue,
label = "Reference", xlabel = "", ylabel = "Precipitation (mm h⁻¹)",
tickfontsize = 12,
guidefontsize = 14, legend_font_pointsize = 12, 
legend_foreground_color = nothing,
xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")))

p2 = plot(Timestamp, Precip_sh.(t), linewidth = 2, color = :red,
label = "Stress", xlabel = "", ylabel = "Precipitation (mm h⁻¹)",
tickfontsize = 12,
guidefontsize = 14, legend_font_pointsize = 12, 
legend_foreground_color = nothing,
xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")))

p3 = plot(sol_df.Timestamp, sol_df.Sᵥ, linewidth = 2, label = "V stage", color = :blue,
xlabel = "", ylabel = "Stage",
tickfontsize = 12,
guidefontsize = 14, legend_font_pointsize = 12, 
legend_foreground_color = nothing,
xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")))
plot!(sol_df.Timestamp, sol_df.Sᵣ, label = "R stage", linewidth = 2)
p3
l = @layout [a{1.0w} ; b{1.0w}; c{1.0w}]


p=plot(p1,p2, p3; layout =l, size = (600, 900), left_margin = (4.0, :mm))

savefig(p, "./img/_Divagro_input.png")

tm_ticks = round.(Timestamp, Month(1)) |> unique
plot(sol_df.Timestamp, sol_df.Sᵥ, linewidth = 2, label = "V stage", color = :blue,
xlabel = "", ylabel = "Stage",
tickfontsize = 12,
guidefontsize = 14, legend_font_pointsize = 12, 
legend_foreground_color = nothing, xlimits = (Timestamp[1], Timestamp[end]),
xticks = (tm_ticks[1:5], Dates.format.(tm_ticks[1:5], "mm-dd")))
plot!(sol_df.Timestamp, sol_df.Sᵣ, label = "R stage", linewidth = 2)
savefig("./img/_Divagro_phenology.png")