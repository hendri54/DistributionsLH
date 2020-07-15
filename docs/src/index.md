# DistributionsLH

This package makes `ModelObject`s for random variables drawn from select distributions.

For the most part, this just wraps distributions from `Distributions.jl` into `ModelObjects` from `ModelParams` so that the parameters can be calibrated conveniently and so that fixed properties can be set using switches objects.

There are also extensions to the distributions in `Distributions.jl`. For example, the `Beta` distribution is defined over an arbitrary range rather than [0, 1].

# Generic interface

Each distribution implements the following

```@docs
init_distribution
draw
scale_normal_draws
quantiles
min_value
max_value
check_draws
```

# Uniform distribution

```@docs
AbstractUniformSwitches
AbstractUniform
Uniform
UniformFixedBounds
UniformCenteredSwitches
UniformCentered
init_uniform
```

# Beta distribution

```@docs
BetaSwitches
Beta
init_beta
```

--------------