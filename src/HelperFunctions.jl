#-------------------------------------------------------------------------#
# Functions to make our life easier and to ensure correct transformations #
#-------------------------------------------------------------------------#

# Unit conversion functions
cm2MPa(x) = 98.1e-6 * x
MPa2cm(x) = 1/98.1e-6 * x

# Smooth minimum function
function smoothmin(x, y; ε = 1e-6)
    return (x + y - sqrt((x - y)^2 + ε)) / 2
end

# Smooth step function: returns 0 near threshold, 1 elsewhere
function smooth_ifelse(x; threshold=-15000, width=500)
    # Normalize x into [0,1] range around threshold
    t = clamp((x - threshold) / width, 0, 1)
    return t^2 * (3 - 2t)  # cubic smoothstep
end

