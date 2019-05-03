using Plots

C(fish, boats) = fish * exp(-47/boats)
effort = 0:100
fish = 0:10:4000
land = [C(f, e) for f in fish, e in effort]

heatmap(effort, fish, land, xlabel="Effort (boats)", ylabel="Biomass (t fish)")

plot(f -> C(f, 100), fish)
plot(e -> C(2000, e), 1:1000)


spillover(R, M, μ, α) = μ * (R - (α / (1-α)) * M)
logistic(n, r, K) = 1 + r * (1 - n/K)
popgrowth(N, r, K, A) = N * logistic(N/A, r, K)
landings(n, boats) = n * exp(-47/boats)


μ = 0.05
r = 0.1
α = 0.1
K = 1.0
A = 100.0
N  = A*K
M = (1-α) * N
R = α * N

plot(N -> popgrowth(N, r, K, A), 0, 1200, ylims=(0, 310))
plot!(N -> popgrowth(N, r, K, A*(1-α)), 0, 1200)
plot!(N -> spillover(R, N, μ, α), 0, 1200)



a = 0.5 # area^2 s^-1
Th = 0.5 # boat days (unit fish)^-1
catch_per_boat(n) = a * n / (1 + a * Th * n)
harvest(n, B) = B * catch_per_boat(n)
plot(n -> catch_per_boat(n), 0, 1000)
plot(n -> harvest(n, 10), 0, 1000)

using Polynomials
B = 3 # boats
roots(Poly([0, r-B*a, a*Th*r-r/K, a*Th*r*K]))

plot(B -> maximum(roots(Poly([0, r-B*a, a*Th*r-r/K, a*Th*r*K]))), 1, 10000)
