include("../src/mpa_model.jl")
m = MpaModel

K = 1 # ton fish km^-2
r = 0.1
A = 100 # km^2
α = 0.
μ = 0.1

# both start at carrying capacity
Nreserve = A * α * K
Nopen = A * (1-α) * K

Nfishers = 20
a = 0.1 # km^2 d^-1
Th = 0.001 # boat d (ton fish)^-1


bfishers = [m.BasicFisher(a, Th) for i in 1:Nfishers]
cfishers = [m.CatchObservingFisher(a, Th, 0.05, rand()*0.5) for i in 1:Nfishers]
mfishers = [m.CatchObservingFisher(a, Th, 0.05, rand()*0.5) for i in 1:Nfishers]
mixedfishers = [bfishers; cfishers; mfishers]

for fishers in [bfishers, cfishers, mfishers, mixedfishers]
    ground = m.FishingGround(K, r, A, α, μ, Nreserve, Nopen)
    m.grow_population!(ground)
    m.spillover!(ground)
    # No growth or spillover when already at carrying capacity?
    @assert ground.Nreserve ≈ Nreserve
    @assert ground.Nopen ≈ Nopen

    m.harvest!(ground, fishers)
    # catch reduces open population but not reserve population
    @assert ground.Nreserve ≈ Nreserve
    @assert 0 <= ground.Nopen < Nopen

    nr1 = ground.Nreserve
    no1 = ground.Nopen
    m.spillover!(ground)
    # spillover transfers fish from reserve to open
    @assert ground.Nreserve <= nr1
    @assert ground.Nopen >= no1
    @assert ground.Nopen - no1 ≈ nr1 - ground.Nreserve


    nr1 = ground.Nreserve
    no1 = ground.Nopen
    ntot = nr1 + no1
    m.set_protected!(ground, ground.α + 0.1)
    @assert ground.Nreserve > nr1
    @assert ground.Nopen <  no1
    @assert ground.Nopen + ground.Nreserve ≈ ntot
    nr1 = ground.Nreserve
    no1 = ground.Nopen
    ntot = nr1 + no1
    m.set_protected!(ground, ground.α - 0.1)
    @assert ground.Nreserve < nr1
    @assert ground.Nopen > no1
    @assert ground.Nopen + ground.Nreserve ≈ ntot


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
end
