using Plots
include("../src/mpa_model.jl")
m = MpaModel

K = 1 # ton fish km^-2
r = 0.1
A = 100 # km^2
α = 0.1
μ = 0.1
# both start at carrying capacity
Nreserve = A * α * K
Nopen = A * (1-α) * K

Nfishers = 20
a = 0.1 # km^2 d^-1
Th = 0.001 # boat d (ton fish)^-1


ground = m.FishingGround(K, r, A, α, μ, Nreserve, Nopen)
fishers = [m.Fisher(a, Th) for i in 1:Nfishers]

m.grow_population!(ground)
m.spillover!(ground)
# No growth or spillover when already at carrying capacity?
@assert ground.Nreserve ≈ Nreserve
@assert ground.Nopen ≈ Nopen

m.catch_fish!(ground, fishers)
# catch reduces open population but not reserve population
@assert ground.Nreserve ≈ Nreserve
@assert 0 <= ground.Nopen < Nopen

nr1 = ground.Nreserve
no1 = ground.Nopen
m.spillover!(ground)
# spillover transfers fish from reserve to open
@assert ground.Nreserve < nr1
@assert ground.Nopen > no1
@assert ground.Nopen - no1 ≈ nr1 - ground.Nreserve


ground.Nreserve = Nreserve
ground.Nopen = Nopen
nsim = 100
NNres = zeros(nsim)
NNopen = zeros(nsim)
landings = zeros(nsim)
for t in 1:nsim
    NNres[t] = ground.Nreserve
    NNopen[t] = ground.Nopen
    landings[t] = m.update!(ground,  fishers)
end

@assert all(landings .<= NNopen)
@assert all(NNres .>= 0)
@assert all(NNopen .>= 0)

plot(NNopen, label="Open biomass")
plot!(NNres, label="MPA biomass")
plot!(landings, label="Landings")

n = 100
DD = range(0, 1, length=50)
TT = range(0, .1, length=50)
C = [m.catch_per_boat(n, D, T) for  D in DD, T in TT]
heatmap(TT, DD, C, xlab="Th", ylab="D")
