using Plots
include("mpa_model.jl")
m = MpaModel

function run_simulation(params, nsim)
    K, r, A, α0, α100, μ, a, Th, αpref, Δα, Nfishers = params
    # Biomass in reserve and open areas. Both start at carrying capacity
    Nreserve = A * α0 * K
    Nopen = A * (1-α0) * K

    # Instantiate an array of OpinionatedFishers
    fishers = [m.OpinionatedFisher(a, Th, αpref, Δα) for i in 1:Nfishers]
    # Instantiate a FishingGround
    ground = m.FishingGround(K, r, A, α0, μ, Nreserve, Nopen)

    NNres = zeros(nsim)
    NNopen = zeros(nsim)
    spill = zeros(nsim)
    landings = zeros(nsim)
    αα_avg = fill(α0, nsim)

    for t in 1:nsim
        NNres[t] = ground.Nreserve
        NNopen[t] = ground.Nopen
        spill[t] = m.spillover(ground.Nreserve, ground.Nopen, ground.μ, ground.α)
        landings[t] = m.update!(ground, fishers)
        if t == 100
            m.set_protected!(ground, α100)
            αα_avg[t] = α100
        end
        if t > 100
            αα_avg[t] = m.α_consensus(fishers)
            m.set_protected!(ground, αα_avg[t])
        end
    end
    return (NNopen, NNres, spill, landings, αα_avg)
end

function plot_results(NNopen, NNres, spill, landings, αα_avg, burnin=1)
    p1 = plot(NNopen[burnin:end], label="Open biomass");
    plot!(p1, NNres[burnin:end], label="MPA biomass");
    plot!(p1, landings[burnin:end], label="Landings");
    plot!(p1, spill[burnin:end], label="Spillover")
    p2 = plot(landings[burnin:end], label="Landings", linecolor=:green);
    p3 = plot(αα_avg[burnin:end], label="MPA fraction", linecolor=:black);
    p = plot(p1, p2, p3, layout=(3,1))
    return p
end

# Fishing ground parameters
K = 1          # carrying capacity, ton fish km^-2
r = 1.5        # Intrinsic growth rate of fish populations
A = 100        # Area of fishing ground, km^2
α0 = 0.0       # Proportion of fishing ground in MPA at start
α100 = 0.3    # Proportion of fishing ground in MPA after 100 timesteps
μ = 0.5        # reserve spillover parameter

# Fisher parameters
a = 10        # Search rate, km^2 d^-1
Th = 1        # Handling time, boat d (ton fish)^-1
αpref = α100    # The fishers' preferred reserve fraction at beginning
Δα = 1          # How fast the fishers adjust their preferred reserve size
Nfishers = 25  # Number of fishers


plot(n -> Nfishers * m.catch_per_boat(n, a, Th), 0, K,
    xlabel="Fish biomass density (T km^-2)", ylabel="Total catch (T)")
plot!(n -> A*n*m.logistic(n, r, K), 0, K)

plot(N -> Nfishers * m.catch_per_boat(N/A, a, Th), 0, A * K,
    xlabel="Total fish biomass (T)", ylabel="Total catch (T)", )
plot!(N -> N * m.logistic(N/A, r, K), 0, A*K)

params1 = (K=K, r=r, A=A, α0=α0, α100=α100, μ=μ,
    a=a, Th=Th, αpref=αpref, Δα=Δα, Nfishers=Nfishers)
nsim = 300
res1 = run_simulation(params1, nsim) #NNopen, NNres, landings, αα_avg
p1 = plot_results(res1..., 1);

params2 = (K=K, r=r, A=A, α0=α0, α100=α100, μ=μ,
    a=a, Th=Th, αpref=0., Δα=Δα, Nfishers=Nfishers)
res2 = run_simulation(params2, nsim)
p2 = plot_results(res2..., 1);
plot(p1, p2)
