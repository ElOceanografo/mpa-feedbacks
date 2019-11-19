module MpaModel
using Statistics

###############################################################################
# Fishing ground
###############################################################################

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

spillover(Nr, No, μ, α) = μ * (Nr - (α / (1-α)) * No)
logistic(n, r, K) = r * (1 - n/K)
catch_per_boat(n, a, Th) = a * n / (1 + a * Th * n)

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

###############################################################################
# Fisher types
###############################################################################
abstract type AbstractFisher end

##############
# Basic fisher
##############
struct BasicFisher{T<:Real} <: AbstractFisher
    a::T
    Th::T
end
BasicFisher(a, Th) = BasicFisher(promote(a, Th)...)
α_opinion(fisher::BasicFisher) = 0
update_opinion!(fisher::BasicFisher) = 0

function harvest!(g::FishingGround, fisher::BasicFisher)
    return catch_per_boat(nopen(g), fisher.a, fisher.Th)
end

#########################
# Catch-observing fisher
#########################
mutable struct CatchObservingFisher{T<:Real} <: AbstractFisher
    a::T
    Th::T
    α::T    # Preferred protected fraction
    Δα::T   # Sensitivity to change in catch
    lastcatch::T # Catch last timestep
    thiscatch::T # Catch this timestep

    function CatchObservingFisher{T}(a::T, Th::T, α::T, Δα::T, lastcatch::T, thiscatch::T) where T<:Real
        0 <= α <= 1 || throw(ArgumentError("α=$α is not between 0 and 1"))
        # 0 <= Δα <= 1 || throw(ArgumentError("Δα=$Δα is not between 0 and 1"))
        return new(a, Th, α, Δα, lastcatch, thiscatch)
    end
end
function CatchObservingFisher(a::T, Th::T, α::T, Δα::T, lastcatch::T, thiscatch::T) where T<:Real
    return CatchObservingFisher{T}(a, Th, α, Δα, lastcatch, thiscatch)
end
CatchObservingFisher(a, Th, α, Δα, lastcatch=0, thiscatch=0) = CatchObservingFisher(promote(a, Th, α, Δα, lastcatch, thiscatch)...)

α_opinion(fisher::CatchObservingFisher) = fisher.α

function harvest!(g::FishingGround, fisher::CatchObservingFisher)
    c = catch_per_boat(nopen(g), fisher.a, fisher.Th)
    fisher.lastcatch = fisher.thiscatch
    fisher.thiscatch = c
    return c
end

function update_opinion!(fisher::CatchObservingFisher)
    fisher.lastcatch <= 0 ? change = 0 : change = (fisher.thiscatch - fisher.lastcatch) / fisher.lastcatch
    fisher.α *= 1+fisher.Δα*change
    fisher.α = max(0, min(fisher.α, 1))
end

##########################
# Reserve-observing fisher
##########################
mutable struct MpaObservingFisher{T<:Real} <: AbstractFisher
    a::T
    Th::T
    α::T    # Preferred protected fraction
    Δα::T   # Sensitivity to change in catch
    lastmpa::T # MPA biomass density last timestep
    thismpa::T # MPA biomass density this timestep

    function MpaObservingFisher{T}(a::T, Th::T, α::T, Δα::T, lastmpa::T, thismpa::T) where T<:Real
        0 <= α <= 1 || throw(ArgumentError("α=$α is not between 0 and 1"))
        # 0 <= Δα <= 1 || throw(ArgumentError("Δα=$Δα is not between 0 and 1"))
        return new(a, Th, α, Δα, lastmpa, thismpa)
    end
end
function MpaObservingFisher(a::T, Th::T, α::T, Δα::T, lastmpa::T, thismpa::T) where T<:Real
    return MpaObservingFisher{T}(a, Th, α, Δα, lastmpa, thismpa)
end
MpaObservingFisher(a, Th, α, Δα, lastmpa=0, thismpa=0) = MpaObservingFisher(promote(a, Th, α, Δα, lastmpa, thicatch)...)

α_opinion(fisher::MpaObservingFisher) = fisher.α

function harvest!(g::FishingGround, fisher::MpaObservingFisher)
    c = catch_per_boat(nopen(g), fisher.a, fisher.Th)
    fisher.lastmpa = fisher.thismpa
    fisher.thismpa = nreserve(g)
    return c
end

function update_opinion!(fisher::MpaObservingFisher)
    fisher.lastmpa <= 0 ? change = 0 : change = (fisher.thismpa - fisher.lastmpa) / fisher.lastmpa
    fisher.α *= 1+fisher.Δα*change
    fisher.α = max(0, min(fisher.α, 1))
end

###############################################################################
# Functions to update whole system
###############################################################################

function harvest!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    n = nopen(g)
    landings = 0
    for f in fishers
        landings += harvest!(g, f)
    end
    landings >= g.Nopen ? g.Nopen = 0 : g.Nopen -= landings
    return landings
end

function α_consensus(fishers::Vector{T}) where T <: AbstractFisher
    return mean(α_opinion.(fishers))
end

function update!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    grow_population!(g)
    spillover!(g)
    landings = harvest!(g, fishers)
    update_opinion!.(fishers)
    set_protected!(g, α_consensus(fishers))
    return landings
end

end # module
