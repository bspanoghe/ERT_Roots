#=--------------------------------------------------------------
Flux matrix potential for integration of hydraulic conductivity 
from bulk soil to soil-root interface.
--------------------------------------------------------------=#
function Khfunc(h, α, n, K_s, l)
    m = 1 - 1 / n
    α_h = α * abs(h)
    q = K_s * (1 - α_h^(n - 1) * (1 + α_h^n)^(-m))^2 *
        (1 + α_h^n)^(-l * m)
    return q
end

function fluxmpfunction(hup, α, n, Ks, l)
    hmin = -100000.0  # cm, lower boundary of integration
    fluxmp = quadgk(h -> Khfunc(h, α, n, Ks, l), hmin, hup)[1]
    return fluxmp
end

function itp_matricflux(h) 
    lfun = linear_interpolation(hrange, MP, extrapolation_bc=Line())
    return lfun(h)
end
#@register_symbolic itp_matricflux(h)

function ksoilfun(hsoil, hint)
    (itp_matricflux(hsoil) - itp_matricflux(hint)) / (hsoil - hint + eps())
end
@register_symbolic ksoilfun(hsoil, hint)


#=------------------------------------------------------------------------------
vanGenuchten(h, θ_s, θ_r, α, n, K_s, l)
vanGenuchten: Van Genuchten Mualem model for unsaturated soil 
van Genuchten, M.Th., 1980. A closed-form equation for predicting the   
hydraulic conductivity of unsaturated soils. Soil Sci. Soc. Am. J. 44,  892898.
-------------------------------------------------------------------------------=#
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


#=---------------------------------------------
Top and bottom boundary conditions functions
----------------------------------------------=#
# Describe the boundary condition at the top (evaporation or precipitation)
# Water fluxes are positive in the upward direction (precipitation is negative and evaporation positive)
function bc_top(htop, Ktop, Ksat, precevap, K, h₁, dz)

    # Determine "rain" (precevap < 0)
    is_rain = precevap < 0  # use broadcasting-safe symbolic comparison

    # htop and Ktop values depending on precipitation or evaporation
    htop_val = ifelse(is_rain, 0.0, htop)
    Ktop_val = ifelse(is_rain, Ksat, Ktop)

    kminus = 0.5 * (K + Ktop_val)
    qtop = -kminus * ((htop_val - h₁) / dz + 1)

    s2 = 0.5*(1 + tanh(10*(h₁ - 1e-6)))
    s1 = 0.5*(1 + tanh(-10*(precevap + 1e-6)))  # smooth transition around zero
    
    # Compute qsup depending on rain/evaporation
    maxinf = qtop * (1 - (s1*s2))
    qsup_precip = ifelse(abs(precevap) > abs(maxinf), maxinf, precevap)
    qsup = ifelse(is_rain, qsup_precip, min(qtop, precevap))

    return qsup
end

#bottom boundary conditions to be done (for now: free darinage)

#=-----------------------------------------
Functions for solving soil-root transport
------------------------------------------=#

B_fun(ρ) = ρ ≈ 0 ? 0.0 : (2* (ρ^2 - 1))/(1 - (0.53*ρ)^2 + 2*ρ^2*(log(0.53) + log(ρ)))
@register_symbolic B_fun(ρ)

rhizo_fun(rld) = rld ≈ 0.0 ? 0.0 : 1/(π * rld)^0.5 
@register_symbolic rhizo_fun(rld)

#=----------------------------------------------------------------
Function to calculate change in water potential using elasticity
----------------------------------------------------------------=#
dh_fun(W, dW, ε, tiny, dh_surr) = W < tiny ? dh_surr : ε/W*dW
@register_symbolic dh_fun(W, dW, ε, tiny, dh_surr)