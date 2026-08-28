function soil_module(; name, α, n, Kₛ, l, θₛ, θᵣ, dz)
    params = @parameters(
        α = α, n = n, Kₛ = Kₛ, l = l, θₛ = θₛ, θᵣ = θᵣ,
        dz = dz,
    )
    vars = @variables (
        h(t)[1:nz], C(t)[1:nz], K(t)[1:nz], θ(t)[1:nz],
        h_top(t), K_top(t),
        prec_evap(t), 
        s(t)[1:nz],
        dh(t)[1:nz],
        K_half(t)[1:(nz)],
        q(t)[1:(nz+1)],
        dq(t)[1:nz], hT(t)[1:nz]
    )
    eqs = [
        [C[i] ~ vanGenuchten_C(h[i], θₛ, θᵣ, α, n) for i in 1:nz]...,
        [K[i] ~ vanGenuchten_K(h[i], θₛ, θᵣ, α, n, Kₛ, l) for i in 1:nz]...,
        [θ[i] ~ vanGenuchten_θ(h[i], θₛ, θᵣ, α, n) for i in 1:nz]...,

        h_top ~ hmin,  # example value for htop
        K_top ~ vanGenuchten_K(h_top, θₛ, θᵣ, α, n, Kₛ, l), # K in

        # Richards equation in 1D
        K_half[1] ~ K_top,
        [K_half[i+1] ~ 2 / (1/K[i] + 1/K[i+1]) for i in 1:(nz-1)]...,
    
        [q[i+1] ~ K_half[i+1] * ( (h[i+1]-h[i]) / dz -  dz/dz ) for i in 1:(nz-1)]...,
    
        # topflux
        q[1] ~ bc_top(h_top, K_top, Kₛ, prec_evap, K[1], h[1], dz),
        # bottomflux free drainage
        q[nz+1] ~ - K[nz], 

        [dq[i] ~ (q[i+1] - q[i])  for i in 1:nz]...,

        # time derivative: C * dh/dt = -dq - s  (signs depend on conventions)
        [dh[i] ~ ( dq[i]/ dz - s[i]/dz ) / C[i] for i in 1:nz]..., 
  
        [D(h[i]) ~ dh[i] for i in 1:nz]...,
        [hT[i] ~ z[i] + h[i] for i in 1:nz]...,        # eq. 7 in paper
    ]

    system = ODESystem(eqs, t; name)
    return system
end

function rootuptake_module(; name,  εₓ, rᵣ, kᵣ, kₓ, Ψ_ref, k_Ψ, kc)
    params = @parameters (εₓ = εₓ, rᵣ = rᵣ, kᵣ = kᵣ, kₓ = kₓ, Ψ_ref = Ψ_ref, k_Ψ = k_Ψ, kc = kc)
    vars = @variables (        
        F(t)[1:nz], Tp(t), Kₓ(t)[1:nz], dz(t), lᵣ(t)[1:nz],
        Hₓ(t)[1:nz], dHₓ(t)[1:nz], dWₓ(t)[1:nz], Wₓ(t)[1:nz],
        hₛ(t)[1:nz], dhₛ(t), Hₛ(t)[1:nz], s(t)[1:nz], r_rhiz(t)[1:nz],
        rld(t)[1:nz], ρ(t)[1:nz], B(t)[1:nz], Hᵣₛ(t)[1:nz],
        Kᵣ(t)[1:nz], uptake(t)[1:nz], H₀(t),
        LAI(t), f_Ψ(t)
    )
    eqs = [
        # root-soil dimension aspects
        [lᵣ[i] ~ rld[i] * dz for i in 1:nz]..., # implicit multiplication with 1 cm² soil surface area
        [r_rhiz[i] ~ rhizo_fun(rld[i]) for i in 1:nz]...,
        [ρ[i] ~ r_rhiz[i] / rᵣ for i in 1:nz]...,
        [B[i] ~ B_fun(ρ[i]) for i in 1:nz]...,

        # Hydraulic conductivities based on root dimensions and intrinsic k's (h⁻¹)
        [Kᵣ[i] ~  kᵣ * (2*π * rᵣ * lᵣ[i]) / rᵣ  for i in 1:nz]...,
        [Kₓ[i] ~ kₓ / dz * (lᵣ[i]/dz * π * rᵣ^2) for i in 1:nz]...,

        # root xylem flows (positive upward) (cm h⁻¹)
        F[1] ~ (Hₓ[1] - H₀)*Kₓ[1],
        [F[i] ~ (Hₓ[i] - Hₓ[i-1])*(Kₓ[i]*Kₓ[i-1]*2/(Kₓ[i] + Kₓ[i-1])) for i in 2:nz]..., 
        
        # solving for root-soil interface water potential (cm)
        [kᵣ * Hₓ[i] + B[i] * Hₛ[i] * ksoilfun(hₛ[i], Hᵣₛ[i] - z[i]) ~
         Hᵣₛ[i] * (B[i] * ksoilfun(hₛ[i], Hᵣₛ[i] - z[i]) + kᵣ)
         for i in 1:nz]...,

        # water uptake s (cm h⁻¹) per soil layer i
        [s[i] ~ (Hᵣₛ[i] - Hₓ[i]) * Kᵣ[i] for i in 1:nz]...,

        # representative "water mass" of the root xylem (g)
        [Wₓ[i] ~ lᵣ[i] * π * ρ_w * rᵣ^2 for i in 1: nz]...,

        # change in water mass of the root xylem (only for change in Hₓ, not for integration!!!)
        dWₓ[nz] ~ s[nz] - F[nz], 
        [dWₓ[i] ~ s[i] - F[i] + F[i+1] for i in 1:(nz-1)]...,

        # change in xylem water potential (cm h⁻¹)
        dHₓ[1] ~ dh_fun(Wₓ[1], dWₓ[1], εₓ, tiny, dhₛ),
        [dHₓ[i] ~ dh_fun(Wₓ[i], dWₓ[i], εₓ, tiny, dHₓ[i-1])  for i in 2:nz]..., 

        # collar water potential (cm)
        H₀ ~ max(Hₓ[1] - Tp/(Kₓ[1] + tiny), hmin),
        #Ψ ~ H₀ * 98.1e-6,
        f_Ψ ~ 0.5 *( 1 + tanh(k_Ψ * (H₀ * 98.1e-6 - Ψ_ref))),
        Tp ~ ET0_inputfun(t) /10.0 * kc * (1 - exp(-0.45*(LAI))) * f_Ψ,
        # differential equations
        [D(Hₓ[i]) ~ dHₓ[i] for i in 1:nz]...,
        [D(uptake[i]) ~ s[i] for i in 1:nz]...,
    ]
    system = ODESystem(eqs, t; name)
    return system
end

function phenology_module(; name, Tmin, Tmax, Topt, v_max, S_ref, k_s, k_Ψ, Ψ_ref, r_LAI, r_max)
    params = @parameters(
        Tmin = Tmin,
        Tmax = Tmax,
        Topt = Topt, 
        v_max = v_max,
        S_ref = S_ref,
        k_s = k_s,
        k_Ψ = k_Ψ,
        Ψ_ref = Ψ_ref,
        r_LAI = r_LAI,
        r_max = r_max,
    )
    vars = @variables (
        f_T(t), 
        f_R(t),
        f_Ψ(t),
        T(t),
        Sᵥ(t),
        dSᵥ(t), 
        LAI(t),
        dLAI(t),
        Ψ(t),
        Sᵣ(t),
        dSᵣ(t),
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