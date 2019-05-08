using Plots
include("mpa_model.jl")
m = MpaModel

K = 1 # ton fish km^-2
r = 0.5
A = 100 # km^2
α = 0.0
μ = 0.01

# both start at carrying capacity
Nreserve = A * α * K
Nopen = A * (1-α) * K

Nfishers = 200
a = 30 # km^2 d^-1
Th = 0.05 # boat d (ton fish)^-1

fishers = [m.BasicFisher(a, Th)]
ground = m.FishingGround(K, r, A, α, μ, Nreserve, Nopen)

nsim = 200
NNres = zeros(nsim)
NNopen = zeros(nsim)
landings = zeros(nsim)
for t in 1:nsim
    NNres[t] = ground.Nreserve
    NNopen[t] = ground.Nopen
    if t == 50
        global α = 0.2
        m.set_protected!(ground, α)
    end
    landings[t] = m.update!(ground, fishers)
end

p1 = plot(NNopen, label="Open biomass");
plot!(p1, NNres, label="MPA biomass");
plot!(p1, landings, label="Landings");

p2 = plot(landings, label="Landings", );
# plot!(p2, NNres, label="MPA biomass");
# hline!(p2, [A * α * K], label="MPA B0")
plot(p1, p2, layout=(2,1))
