#=---------------------------------------------------------------------------
Root growth functions to generate an rld profile in z and t (own equations)
---------------------------------------------------------------------------=#
# root depth evolution as a function of t (Yin function)
function rootgrowth(t, tₘ, tₑ, dₑ, sowdepth)
    t= smoothmin(t, tₑ)
    root_depth = dₑ*(1 + (tₑ - t)/(tₑ - tₘ))*(t/tₑ)^(tₑ/(tₑ - tₘ)) + sowdepth
    return root_depth
end

# normalized rld profile - p_r is relative position: 0 is root tip, 1 is sowing depth
function rld_norm(p_r, p_max, p_half)
    p_r = smoothmin(1, p_r)
    rld_prof = 1.0 * (1 + ((1-p_max) - (1-p_r))/((1-p_max) - (1-p_half)))*((1-p_r)/(1-p_max))^((1-p_max)/((1-p_max) - (1-p_half)))
    return rld_prof
end

# the maximum rld proportional to root depth d
function rld_max(d, rld_d)
    rld_max = d * rld_d
    return rld_max
end

# rld profile in depth (z) and time (t)
function root_fun(z, t, tₘ, tₑ, dₑ, sow_depth, rld_d, p_max, p_half)
    root_depth = rootgrowth(t, tₘ, tₑ, dₑ, sow_depth)
    rld_m = rld_max(root_depth, rld_d)
    rld = rld_norm(z/root_depth, p_max, p_half)*rld_m
    return rld
end