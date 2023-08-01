# DistributionsLH

Distributions with calibrated parameters.

The main purpose is to wrap distributions into containers so that their parameters can be calibrated using the interface defined in `ModelParams`. Distributions are `ModelObject`s.

[Documentation](lhendricks.org/julia/DistributionsLH/index.html).

## Change Log 2023

Aug-1: ModelParams v4
May-4: `init_normal_switches` accepts symbols for mean and std. Avoids unicode characters that cause issues with latex.
Apr-5: Added DegenerateDistribution.
Feb-3: Replaced ModelParams.value with pvalue.

# ----------------
