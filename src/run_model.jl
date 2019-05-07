using Plots
include("../src/mpa_model.jl")
m = MpaModel

K = 1 # ton fish km^-2
r = 0.1
A = 100 # km^2
α = 0.0
μ = 0.2

# both start at carrying capacity
Nreserve = A * α * K
Nopen = A * (1-α) * K

Nfishers = 100
a = 10 # km^2 d^-1
Th = 0.001 # boat d (ton fish)^-1

fishers = [m.OpinionatedFisher(a, Th, 0.05, rand()*0.5)]
ground = m.FishingGround(K, r, A, α, μ, Nreserve, Nopen)

nsim = 100
NNres = zeros(nsim)
NNopen = zeros(nsim)
landings = zeros(nsim)
for t in 1:nsim
    NNres[t] = ground.Nreserve
    NNopen[t] = ground.Nopen
    if t == 50 m.set_protected!(ground, 0.2) end
    landings[t] = m.update!(ground, fishers)
end

@assert all(landings .<= NNopen)
@assert all(NNres .>= 0)
@assert all(NNopen .>= 0)

p1 = plot(NNopen, label="Open biomass")
plot!(p1, NNres, label="MPA biomass")
plot!(p1, landings, label="Landings")
p2 = plot(landings, label="Landings")
plot(p1, p2, layout=(2,1))
