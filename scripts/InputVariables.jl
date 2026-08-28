#=-----------------------------------------------
Read the input data shared by Jan Vanderborght 
This script is called by main.jl    
Requires InputParameters.jl to be run first           
#-----------------------------------------------=#

# ET0 - reference evapotranspiration (cm/day) at hourly intervals
ET0_df = CSV.read("./data/ET0.dat", DataFrame; delim = ' ', header = true, ignorerepeated=true)
rename!(ET0_df, [:Time, :ET0])

# Precipitation data (cm/day) at daily intervals
Precip_df = CSV.read("./data/Precip.dat", DataFrame; delim = ' ', header = true, ignorerepeated=true)
rename!(Precip_df, [:Time, :Precip])

# Biomass data (g/plant) and crop coefficent (Kc) at weekly intervals - used for transpiration calculation
Biomass_df = CSV.read("./data/BiomassdataLM_logistic.dat", DataFrame; delim = ' ', 
                      header = false, skipto = 2, ignorerepeated=true)
rename!(Biomass_df, [:Time, :Biomass, :kc])

#----------------------------------------------------------------------#
# Combine and process the input data regarding units and time indexing #
#----------------------------------------------------------------------#
# join ET0 (daily) with precipitation which is indexed by Time
transform!(ET0_df, :Time => ByRow(t -> floor.(t)) => :Day)
ET0_df = leftjoin(ET0_df, Precip_df, on = Pair(:Day, :Time))
ET0_df = subset(ET0_df, :Time => t -> t .< 157)

# Convert precipitation and ET0 from cm/day to cm/hour
transform!(ET0_df, :Precip => ByRow(x -> x/(24.0)) => :Precip) 
transform!(ET0_df, :ET0 => ByRow(x -> x/(24.0)) => :ET0) 
transform!(ET0_df, :Time => ByRow(t -> t*24.0) => :Time)  # convert time from days to hours

# Convert Biomass time from weeks to hours
transform!(Biomass_df, :Time => ByRow(t -> t*7.0*24.0) => :Time)

#------------------------------------------------------------------------------#
# Define input functions for ET0, Precipitation, Transpiration and Evaporation #
#------------------------------------------------------------------------------#
function ET0_inputfun(t)
    lfun = LinearInterpolation(ET0_df.Time, ET0_df.ET0)
    return lfun(t)
end

function Biomass_inputfun(t)
    lfun = LinearInterpolation([0, Biomass_df.Time...], [0, Biomass_df.Biomass...])
    return lfun(t)
end

function Precip_inputfun(t)
    lfun = linear_interpolation(ET0_df.Time, ET0_df.Precip)
    return lfun(t)
end

LAI(t) = -1.3126E-5 * (Biomass_inputfun(t) * Plant_density)^2 + 1.1517E-2 * Biomass_inputfun(t) *Plant_density
factor_T(t) = 1-exp(-0.45*(LAI(t)))

Transpiration(t) = ET0_inputfun(t) * kc * factor_T(t)
Evaporation(t) = ET0_inputfun(t) * kc * (1 - factor_T(t))


#---------------------------------------------------------------------------#
# Register the required input functions for the model as symbolic functions #
#---------------------------------------------------------------------------#
@register_symbolic Precip_inputfun(t)
@register_symbolic Transpiration(t)
@register_symbolic Evaporation(t)


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