#---------------------------------------------#
# Visualize the input functions defined above #
#---------------------------------------------#
timerange = collect(800:1:4200)

plot(timerange, ET0_inputfun.(timerange) )
plot(timerange, Precip_inputfun.(timerange))

plot(timerange, Transpiration.(timerange))
plot!(timerange, Evaporation.(timerange), xlims = (0, 2400))

plot(timerange, Biomass_inputfun.(timerange))
plot(timerange, LAI.(timerange))

 tanh(a_g * (P - P_ref)) 

 P_ref = -1.0
 P = collect(-1.70:0.1:0)
 a_g = 3.0
plot(P, 0.5 .+ 0.5 .*tanh.(a_g * (P .- P_ref)))
plot(P, 0.5 .+ 0.5 .*tanh.(a_g * (P_ref .- P)))

log(2) / (log(40- 5/25 - 5))
T_min = 5
T_opt = 25
T_max = 40

α = log(2) / (log(T_max - T_min/T_opt - T_min))


T = collect(0.1:0.1:50)

T_eff = 
TU = (((T_max .- T)/(T_max-T_opt)).*((T .- T_min)/(T_opt-T_min)).^((T_opt-T_min)/(T_max-T_opt))).^1.0

plot(T, TU)

Ψ_ref = -0.75
k_Ψ = 4.0
ψ = collect(-1.5:0.1:0)
plot(ψ, 0.5 .*( 1 .+ tanh.(k_Ψ .* (ψ .- Ψ_ref))))
plot!(ψ, 0.5 .*( 1 .+ tanh.(5.0 .* (ψ .- -1.0))))
