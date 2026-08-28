#----------------------------------------#
# Fixed inputs for the model simulations #
#----------------------------------------#

# mathematical convenience factors 
tiny = 1e-6
hmin = -15000 + 1
hmax = -5 - 1

# soil parameters, discretization and initial conditions
finesoil = (α = 0.0083, n = 1.2539, K_s = 2.272 / (24), l = 0.5, θ_s = 0.43, θ_r = 0.078)

nz = 10  # number of soil layers
Δz = 7 # cm, thickness of each soil layer
depth_hini = [-75, -60, -40, -20, 0]
hini = [-43, -43, -92, -168, -168]


# root system parameters
r_root = 0.05 # cm, root radius 
εₓ = MPa2cm(10.0) # cm (calculated from MPa), xylem elatsic modulus
k_r = 0.1 # cm⁻¹ h⁻¹, radial intrinsic root hydraulic conductivity 
k_x = 0.1 # cm⁻¹ h⁻¹, axial intrinisc root hydraulic conductivity 

t_m = 1200 # h (shape factor of root depth development)
t_e = 1800 # h (time of maximum root depth)
d_e = 65 # cm (maximum rooting depth)
d_s = 2 # cm (sowing depth)
rld_d = 0.02 # (proportionality factor of rld_max to root depth d)
p_max = 0.35 # relative position of maximum rld in rld distribution (0= sow_depth, 1 = root depth)
p_half = 0.7 # shape factor of rld distribution

# Settings for transpiration and evaporation calculation
Plant_density = 9.5 # plants per m2
kc = 1.0 # crop coefficient for calculating potential transpiration


