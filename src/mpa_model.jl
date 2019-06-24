module MpaModel
using Statistics
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
struct BasicFisher{T<:Real} <: AbstractFisher
    D::T
    Th::T
end
BasicFisher(D, Th) = BasicFisher(promote(D, Th)...)
α_opinion(fisher::BasicFisher) = 0
update_opinion!(fisher::BasicFisher, c::Real) = 0

mutable struct OpinionatedFisher{T<:Real} <: AbstractFisher
    D::T
    Th::T
    α::T
    Δα::T
    last::T

    function OpinionatedFisher{T}(D::T, Th::T, α::T, Δα::T, last::T) where T<:Real
        0 <= α <= 1 || throw(ArgumentError("α=$α is not between 0 and 1"))
        # 0 <= Δα <= 1 || throw(ArgumentError("Δα=$Δα is not between 0 and 1"))
        return new(D, Th, α, Δα, last)
    end
end
function OpinionatedFisher(D::T, Th::T, α::T, Δα::T, last::T) where T<:Real
    return OpinionatedFisher{T}(D, Th, α, Δα, last)
end
OpinionatedFisher(D, Th, α, Δα, last=0) = OpinionatedFisher(promote(D, Th, α, Δα, last)...)

α_opinion(fisher::OpinionatedFisher) = fisher.α

function update_opinion!(fisher::OpinionatedFisher, c::Real)
    fisher.last <= 0 ? change = 0 : change = (c - fisher.last) / fisher.last
    fisher.α *= 1+fisher.Δα*change
    fisher.α = max(0, min(fisher.α, 1))
    # if c > fisher.last
    #     fisher.α = min(1, fisher.α + fisher.Δα)
    # else
    #     fisher.α = max(0, fisher.α - fisher.Δα)
    # end
    fisher.last = c
    return fisher.α
end

function α_consensus(fishers::Vector{T}) where T <: AbstractFisher
    return mean(α_opinion.(fishers))
end

spillover(R, M, μ, α) = μ * (R - (α / (1-α)) * M)
logistic(n, r, K) = r * (1 - n/K)
catch_per_boat(n, D, Th) = D * n / (1 + D * Th * n)

# Fish density in reserve/open area. Checks to avoid division by 0
nreserve(g::FishingGround) = g.α > 0 ? g.Nreserve / (g.A * g.α) : 0
nopen(g::FishingGround) = g.α < 1 ? g.Nopen / (g.A * (1-g.α)) : 0

function grow_population!(g::FishingGround)
    nres = nreserve(g)
    nop = nopen(g)
    g.Nreserve += nres * logistic(nres, g.r, g.K) * g.A * g.α
    g.Nopen += nop  * logistic(nop, g.r, g.K) * g.A * (1-g.α)
    return g
end

function spillover!(g::FishingGround)
    s = spillover(g.Nreserve, g.Nopen, g.μ, g.α)
    g.Nreserve -= s
    g.Nopen += s
    return s
end

function catch_fish!(fisher::AbstractFisher, n)
    c = catch_per_boat(n, fisher.D, fisher.Th)
    update_opinion!(fisher, c)
    return c
end

function harvest!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    n = g.Nopen / (g.A * (1-g.α))
    landings = [catch_fish!(f, n) for f in fishers]
    total = sum(landings)
    total >= g.Nopen ? g.Nopen = 0 : g.Nopen -= total
    return total
end

function set_protected!(g::FishingGround, α::Real)
    0 < α < 1 || DomainError(α, "α must be between 0 and 1")
    Δarea = g.A * α - g.A * g.α
    if Δarea > 0
        ΔN = nopen(g) * Δarea
    else
        ΔN = nreserve(g) * Δarea
    end
    g.Nreserve += ΔN
    g.Nopen -= ΔN
    g.α = α
    return g
end

function update!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    grow_population!(g)
    spillover!(g)
    landings = harvest!(g, fishers)
    return landings
end

end # module
