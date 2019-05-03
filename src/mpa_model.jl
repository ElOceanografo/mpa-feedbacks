
module MpaModel
# using DifferentialEquations

mutable struct FishingGround{T<:Real}
    K::T
    r::T
    A::T
    α::T
    μ::T
    Nreserve::T
    Nopen::T
end
FishingGround(K, r, A, α, μ, Nreserve, Nopen) = FishingGround(promote(K, r, A, α, μ, Nreserve, Nopen)...)

abstract type AbstractFisher end
struct Fisher{T<:Real} <: AbstractFisher
    D::T
    Th::T
end
Fisher(D, Th) = Fisher(promote(D, Th)...)

spillover(R, M, μ, α) = μ * (R - (α / (1-α)) * M)
logistic(n, r, K) = 1 + r * (1 - n/K)
popgrowth(N, r, K, A) = N * logistic(N/A, r, K)
catch_per_boat(n, D, Th) = D * n / (1 + D * Th * n)

function grow_population!(g::FishingGround)
    g.Nreserve = logistic(g.Nreserve/(g.A * g.α), g.r, g.K) * g.Nreserve
    g.Nopen = logistic(g.Nopen/(g.A * (1-g.α)), g.r, g.K) * g.Nopen
    return g
end

function spillover!(g::FishingGround)
    s = spillover(g.Nreserve, g.Nopen, g.μ, g.α)
    g.Nreserve -= s
    g.Nopen += s
    return g
end

function catch_fish!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    n = g.Nopen / (g.A * (1-g.α))
    landings = [catch_per_boat(n, f.D, f.Th) for f in fishers]
    total = sum(landings)
    if total >= g.Nopen
        g.Nopen = 0
    else
        g.Nopen -= total
    end
    return total
end

function update!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    grow_population!(g)
    spillover!(g)
    landings = catch_fish!(g, fishers)
    return landings
end

function set_protected!(g::FishingGround, α::Real)
    0 <= α <= 1 || DomainError(α, "α must be between 0 and 1")
    g.α = α
end

end # module
