module DistributionsLH

using DocStringExtensions, Lazy, Parameters, Random, StatsBase, Statistics
using Distributions, LinearAlgebra
using ModelObjectsLH, ModelParams

export AbstractDistributionSwitches, AbstractDistributionLH

# Generic functions
export init_distribution, draw, check_draws, isbounded, min_value, max_value, quantiles, scale_normal_draws
# Uniform
export AbstractUniform, AbstractUniformSwitches, Uniform, UniformSwitches, UniformCentered, UniformCenteredSwitches, UniformFixedBounds
export init_uniform
# Beta
export AbstractBeta, AbstractBetaSwitches, Beta, BetaSwitches
export init_beta, alpha_value, beta_value
# Normal
export NormalSwitches, Normal
export init_normal

abstract type AbstractDistributionSwitches{T1} <: ModelSwitches end
abstract type AbstractDistributionLH{T1} <: ModelObject end

include("uniform.jl")
include("beta.jl");
include("normal.jl");
include("mv_normal.jl")


## --------  Generic

"""
	init_distribution(objId, switches)

Initialize a distribution from its switches.
"""
function init_distribution end

"""
    draw(distrib, nDims, rng)
    draw(distrib, rng)

Draw from the distribution.
"""
function draw end

"""
    check_draws(distrib, drawM)

Check draws against bounds (if any) and `Inf` values.
"""
function check_draws end

"""
	quantiles(distrib, pctV)

Quantiles.
"""
function quantiles end


"""
	$(SIGNATURES)

Converts N(0,1) draws into the specified distribution.
"""
function scale_normal_draws(b :: AbstractDistributionLH{F},  inM :: T)  where {F, T}
    cdfM = cdf.(Distributions.Normal(), inM);
    drawM = quantiles(b,  cdfM);
    return drawM
end


"""
    isbounded(distrib)

Is this a bounded distribution?
"""
function isbounded end

"""
    min_value(distrib)

Minimum value (if any).
"""
function min_value end

"""
    max_value(distrib)

Maximum value (if any).
"""
function max_value end


function pretty_print(x)
    show(stdout, MIME("text/plain"), x);
    println("");
end


end # module
