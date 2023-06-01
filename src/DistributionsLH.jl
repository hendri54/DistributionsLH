module DistributionsLH

using DocStringExtensions, Lazy, Parameters, Random, StatsBase, Statistics
using Distributions, LinearAlgebra
using ModelObjectsLH, ModelParams

export AbstractDistributionSwitches, AbstractDistributionLH

# Generic functions
export init_distribution, draw, check_draws, isbounded, min_value, max_value, quantiles, scale_normal_draws
# Degenerate
export DegenerateSwitches, DegenerateDistribution;
# Uniform
export AbstractUniform, AbstractUniformSwitches, Uniform, UniformSwitches, UniformCentered, UniformCenteredSwitches, UniformFixedBounds
export init_uniform
# Beta
export AbstractBeta, AbstractBetaSwitches, Beta, BetaSwitches
export init_beta, alpha_value, beta_value
# Normal
export NormalSwitches, Normal, LogNormalSwitches, LogNormal
export init_normal, init_log_normal, init_normal_switches, init_log_normal_switches

abstract type AbstractDistributionSwitches{T1} <: ModelSwitches end
abstract type AbstractDistributionLH{T1} <: ModelObject end

include("degenerate.jl");
include("uniform.jl");
include("beta.jl");
include("normal.jl");
include("log_normal.jl");
include("mv_normal.jl")


## --------  Generic

# Makes it easy to implement StructLH.describe
describe(switches :: AbstractDistributionSwitches{T1}) where T1 = 
    ["Distribution"  "$(typeof(switches))"];

function calibrated_string(val, isCalibrated :: Bool)
    if isCalibrated
        return "calibrated";
    else
        return "fixed at $val";
    end
end

function calibrated_string(p :: Param)
    calibrated_string(default_value(p), is_calibrated(p));
end

function calibrated_string(pv :: ParamVector, vName :: Symbol)
    calibrated_string(retrieve(pv, vName));
end



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
