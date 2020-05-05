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
    subsidy::T
    subsidy_α::T
end

function FishingGround(K, r, A, α, μ, Nreserve, Nopen, subsidy=0, subsidy_α=0)
    return FishingGround(promote(K, r, A, α, μ, Nreserve, Nopen, subsidy, subsidy_α)...)
end

spillover(Nr, No, μ, α) = μ * (Nr - (α / (1-α)) * No)
logistic(n, r, K) = r * (1 - n/K)
catch_per_boat(n, a, Th) = a * n / (1 + a * Th * n)

# Fish density in reserve/open area. Checks to avoid division by 0
nreserve(g::FishingGround) = g.α > 0 ? g.Nreserve / (g.A * g.α) : 0
nopen(g::FishingGround) = g.α < 1 ? g.Nopen / (g.A * (1-g.α)) : 0
growthopen(g::FishingGround) = nopen(g) * logistic(nopen(g), g.r, g.K) * g.A * (1-g.α)
growthreserve(g::FishingGround) = nreserve(g) * logistic(nreserve(g), g.r, g.K) * g.A * g.α

function grow_population!(g::FishingGround)
    Δres = growthreserve(g)
    Δop = growthopen(g)
    g.Nreserve = g.Nreserve + Δres
    g.Nopen = g.Nopen + Δop
    return (Δres, Δop)
end

spillover(g::FishingGround) = spillover(g.Nreserve, g.Nopen, g.μ, g.α)
function spillover!(g::FishingGround)
    s = spillover(g)
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
    return α
end

function subsidy(g::FishingGround{T}) where T
    if g.α > g.subsidy_α
        return g.subsidy
    else
        return zero(T)
    end
end


###############################################################################
# Fisher types
###############################################################################
abstract type AbstractFisher end

update_observation!(fisher::T, g::FishingGround) where T<:AbstractFisher = zero(T)
update_catch!(fisher::T, c) where T<:AbstractFisher = zero(T)
update_opinion!(fisher::T) where T<: AbstractFisher = zero(T)
α_opinion(fisher::T) where T<: AbstractFisher = zero(T)

##############
# Basic fisher
##############
struct BasicFisher{T<:Real} <: AbstractFisher
    a::T
    Th::T
end
BasicFisher(a, Th) = BasicFisher(promote(a, Th)...)

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

function update_catch!(fisher::CatchObservingFisher, c)
    fisher.lastcatch = fisher.thiscatch
    fisher.thiscatch = c
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
MpaObservingFisher(a, Th, α, Δα, lastmpa=zero(a), thismpa=zero(a)) = MpaObservingFisher(promote(a, Th, α, Δα, lastmpa, thismpa)...)

α_opinion(fisher::MpaObservingFisher) = fisher.α

function update_observation!(fisher::MpaObservingFisher, g::FishingGround)
    fisher.lastmpa = fisher.thismpa
    fisher.thismpa = nreserve(g)
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
    catches = [catch_per_boat(nopen(g), f.a, f.Th) for f in fishers]
    landings = sum(catches)
    if landings > g.Nopen
        catches .-= (landings - g.Nopen) / length(catches)
        landings = g.Nopen
    end
    g.Nopen -= landings
    for (f, c) in zip(fishers, catches)
        update_catch!(f, c + subsidy(g))
    end
    return landings
end

function observe!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    for f in fishers
        update_observation!(f, g)
    end
end

function α_consensus(fishers::Vector{T}) where T <: AbstractFisher
    return mean(α_opinion.(fishers))
end

function update!(g::FishingGround, fishers::Vector{T}) where T <: AbstractFisher
    grow_population!(g)
    spillover!(g)
    landings = harvest!(g, fishers)
    observe!(g, fishers)
    update_opinion!.(fishers)
    set_protected!(g, α_consensus(fishers))
    return landings
end

end # module
