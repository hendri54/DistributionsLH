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

# Multivariate Normal Distribution

The main purpose is to be able to draw random variables from the "weight matrix" which is basically a lower triangular decomposition of the covariance matrix, but with ones on the diagonal. The draws are then scaled to match the target means and standard deviations of the marginals.

Also computes conditional distributions.

```@docs
MvNormalLHSwitches
MvNormalLH
cov_matrix
check_cov_matrix
draw_from_weights
check_weight_matrix
conditional_distrib
cond_mean_weights
```

--------------