module DistributionsLH

using DocStringExtensions, Parameters, Random
using Distributions
using ModelParams

export AbstractDistributionSwitches, AbstractDistributionLH

# Generic functions
export draw, check_draws, min_value, max_value, quantiles
# Uniform
export AbstractUniform, AbstractUniformSwitches, Uniform, UniformSwitches, UniformCentered, UniformCenteredSwitches, UniformFixedBounds
export init_uniform

abstract type AbstractDistributionSwitches end
abstract type AbstractDistributionLH <: ModelObject end

include("uniform.jl")

end # module
