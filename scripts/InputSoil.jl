hmin = -15000 + 1
hmax = -5 - 1

nz = 10  # number of soil layers
Δz = 7 # cm, thickness of each soil layer
depth_hini = [-75, -60, -40, -20, 0]
hini = [-43, -43, -92, -168, -168] # initial pressure head [cm] at the depth_hini

# Automatic calculations: not to be changed!!!
z = collect(-Δz/2:-Δz:(-nz*Δz))  # elevation head at center of each layer
itp_h0 = linear_interpolation(depth_hini, hini)
h0 = itp_h0.(z)

hrange = hmin:10:hmax
MP = fluxmpfunction.(hrange, soilpar.α, soilpar.n, soilpar.K_s, soilpar.l)
