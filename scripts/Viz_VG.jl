pars_finesoil = (α = 0.0083, n = 1.2539, K_s = 2.272 / (24), l = 0.5, θ_s = 0.43, θ_r = 0.078)
pars_rockwool = (α = 0.1133, n = 4.203, K_s = 450.0, l = 0.5, θ_s = 0.9, θ_r = 0.0145)
pars_sandyloam = (θs = 0.42, θr = 0.07, α  = 0.07, n  = 1.8, Ks = 150 /(24), l  = 0.5)


h_s = -15000:1:-0.001
θ_finesoil = vanGenuchten_θ.(h_s, pars_finesoil.θ_s, pars_finesoil.θ_r, pars_finesoil.α, pars_finesoil.n)
θ_rockwool = vanGenuchten_θ.(h_s, pars_rockwool.θ_s, pars_rockwool.θ_r, pars_rockwool.α, pars_rockwool.n)
θ_sandyloam = vanGenuchten_θ.(h_s, pars_sandyloam.θs, pars_sandyloam.θr, pars_sandyloam.α, pars_sandyloam.n)

plot(θ_finesoil, -h_s ./ 10000,  yscale = :log10, label = "Fine soil", 
ylabel = "Soil water potential (-MPa)", 
xlabel = "Volumetric moisture content (-)", 
legend_foreground_color = nothing, size = (600, 400), dpi=300,
linewidth = 3, tickfontsize = 12, legend_font_pointsize = 12, guidefontsize = 14)
plot!(θ_sandyloam, -h_s ./ 10000,  yscale = :log10, 
label = "Sandy loam", ylabel = "Soil water potential (-MPa)", 
xlabel = "Volumetric moisture content (-)", linewidth = 3)
plot!(θ_rockwool, -h_s ./ 10000,  yscale = :log10, 
label = "Rockwool", ylabel = "Soil water potential (-MPa)", 
xlabel = "Volumetric moisture content (-)", linewidth = 3)
savefig("./img/VG_θ_logscale.png")

plot(θ_finesoil, -h_s ./ 10000, label = "Fine soil", 
ylabel = "Soil water potential (-MPa)", 
xlabel = "Volumetric moisture content (-)", 
legend_foreground_color = nothing, size = (600, 400), dpi=300,
linewidth = 3, tickfontsize = 12, legend_font_pointsize = 12, guidefontsize = 14)
plot!(θ_sandyloam, -h_s ./ 10000,  
label = "Sandy loam", ylabel = "Soil water potential (-MPa)", 
xlabel = "Volumetric moisture content (-)", linewidth = 3)
plot!(θ_rockwool, -h_s ./ 10000,  
label = "Rockwool", ylabel = "Soil water potential (-MPa)", 
xlabel = "Volumetric moisture content (-)", linewidth = 3)
savefig("./img/VG_θ.png")


plot(θ_finesoil, -h_s, yscale = :log10, label = "Fine soil", xlabel = "Pressure head (cm)", ylabel = "Volumetric moisture content (-)")
plot!(θ_rockwool, -h_s, yscale = :log10, label = "Rockwool", xlabel = "Pressure head (cm)", ylabel = "Volumetric moisture content (-)")

K_finesoil = vanGenuchten_K.(h_s, pars_finesoil.θ_s, pars_finesoil.θ_r, pars_finesoil.α, pars_finesoil.n, pars_finesoil.K_s, pars_finesoil.l)
K_rockwool = vanGenuchten_K.(h_s, pars_rockwool.θ_s, pars_rockwool.θ_r, pars_rockwool.α, pars_rockwool.n, pars_rockwool.K_s, pars_rockwool.l)    
K_sandyloam = vanGenuchten_K.(h_s, pars_sandyloam.θs, pars_sandyloam.θr, pars_sandyloam.α, pars_sandyloam.n, pars_sandyloam.Ks, pars_sandyloam.l)

plot(θ_finesoil, K_finesoil, label = "Fine soil", 
xlabel = "Volumetric moisture content (-)", ylabel = "Hydraulic conductivity (cm/h)",
legend_foreground_color = nothing, size = (600, 400), dpi=300,
linewidth = 3, tickfontsize = 12, legend_font_pointsize = 12, guidefontsize = 14)
savefig("./img/VG_K_clay.png")
plot!(θ_sandyloam, K_sandyloam, label = "Sandy loam", 
xlabel = "Volumetric moisture content (-)", ylabel = "Hydraulic conductivity (cm/h)", linewidth = 3)
savefig("./img/VG_K_2soils.png")

plot!(θ_rockwool, K_rockwool,  label = "Rockwool", 
xlabel = "Volumetric moisture content (-)", ylabel = "Hydraulic conductivity (cm/h)",
 linewidth = 3)

savefig("./img/VG_K_θ.png")

plot(-h_s, K_finesoil, xscale = :log10, label = "Fine soil", xlabel = "Pressure head (cm)", ylabel = "Hydraulic conductivity (cm/h)")
plot!(-h_s, K_rockwool, xscale = :log10, label = "Rockwool", xlabel = "Pressure head (cm)", ylabel = "Hydraulic conductivity (cm/h)")   