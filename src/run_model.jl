using Plots
using Optim

include("mpa_model.jl")
m = MpaModel

# Fishing ground parameters
K = 1          # carrying capacity, ton fish km^-2
r = 1.5        # Intrinsic growth rate of fish populations
A = 100        # Area of fishing ground, km^2
α0 = 0.01       # Proportion of fishing ground in MPA at start
α1 = 0.3    # Proportion of fishing ground in MPA after 100 timesteps
μ = 0.1        # reserve spillover parameter

# Fisher parameters
a = 10        # Search rate, km^2 d^-1
Th = 1        # Handling time, boat d (ton fish)^-1
αpref = α1    # The fishers' preferred reserve fraction at beginning
Δα = 1          # How fast the fishers adjust their preferred reserve size
Nfishers = 25  # Number of fishers

###############################################################################
# Simulations
###############################################################################
# Biomass in reserve and open areas. Both start at carrying capacity
Nreserve = A * α0 * K
Nopen = A * (1-α0) * K


# Instantiate an array of OpinionatedFishers
fishers = [m.BasicFisher(a, Th) for i in 1:Nfishers]
# Instantiate a FishingGround
ground = m.FishingGround(K, r, A, α0, μ, Nreserve, Nopen, 1, 0.05)

m.spillover!(ground)
m.grow_population!(ground)

sum(m.catch_per_boat(m.nopen(ground), f.a, f.Th) for f in fishers)
landings = m.harvest!(ground, fishers)

ground.Nopen
ground.Nreserve

m.subsidy(ground)
