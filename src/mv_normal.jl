export MvNormalLH, MvNormalLHSwitches
export means, stds, n_vars
export cov_matrix, check_cov_matrix, draw_from_weights, check_weight_matrix
export conditional_distrib, cond_mean_weights

"""
	$(SIGNATURES)

Switches used to construct Multivariate Normal object.
"""
mutable struct MvNormalLHSwitches{T1} <: AbstractDistributionSwitches{T1}
    meanV :: Vector{T1}
    stdV :: Vector{T1}
end

"""
	$(SIGNATURES)

Multivariation Normal object. 

This is a `ModelObject` but currently does not support calibrating parameters.
"""
mutable struct MvNormalLH{T1} <: AbstractDistributionLH{T1}
    switches :: MvNormalLHSwitches{T1}
end

MvNormalLH(meanV :: AbstractVector{T1}, stdV :: AbstractVector{T1}) where T1 =
    MvNormalLH(MvNormalLHSwitches(meanV, stdV));

Base.show(io :: IO, switches :: MvNormalLHSwitches{T1}) where T1 =
    print(io, typeof(switches), " with means $(switches.meanV) and stds $(switches.stdV)");

Base.show(io :: IO, m :: MvNormalLH{T1}) where T1 = 
    print(io, typeof(m), " with means $(means(m)) and stds $(stds(m))");

means(m :: MvNormalLH{T1}) where T1 = m.switches.meanV;
stds(m :: MvNormalLH{T1}) where T1 = m.switches.stdV;
n_vars(m :: MvNormalLH{T1}) where T1 = length(means(m));
ModelParams.has_pvector(m :: MvNormalLH{T1}) where T1 = false;

# Make a Distributions.MvNormal that one can sample from
function make_mv_normal(m :: MvNormalLH{T1}, wtM :: AbstractMatrix{T1}) where T1
    return MvNormal(means(m), cov_matrix(m, wtM));
end


"""
	$(SIGNATURES)

Make cov matrix from weight matrix.

The weight matrix is lower triangular with ones on the diagonal.
"""
function cov_matrix(m :: MvNormalLH{T1}, wtM :: AbstractMatrix{T1}) where T1
    @assert check_weight_matrix(m, wtM)
    stdV = stds(m);

    # Scale the weight matrix so that it is the Cholesky decomposition of the Cov matrix
    wt2M = similar(wtM);
    for i1 = 1 : n_vars(m)
       wt2M[i1, :] = wtM[i1, :] ./ sqrt(sum(wtM[i1, :] .^ 2)) .* stdV[i1];
    end
    
    covMatM = wt2M * (wt2M');
    @assert check_cov_matrix(m, covMatM);
    return covMatM
end


"""
	$(SIGNATURES)

Check that weight matrix is valid.
"""
function check_weight_matrix(m :: MvNormalLH{T1}, wtM :: AbstractMatrix{T1}) where T1
    stdV = stds(m);
    n = length(stdV);

    isValid = true;
    isequal(size(wtM), (n, n))  ||  (isValid = false);
    # Lower triangular?
    isValid = isValid  &&  istril(wtM);
    # Diagonal == 1?
    isValid = isValid  &&  all(isapprox.(diag(wtM), one(T1)));
    # If std deviation is 0, variable must have 0 weights
    for i1 = 1 : n
        if stdV[i1] == zero(T1)
            idxV = trues(n);
            idxV[i1] = false;
            isValid = isValid  &&  all(wtM[i1, idxV] .== zero(T1))  &&  
                all(wtM[idxV, i1] .== zero(T1));
        end
    end
    return isValid
end


"""
	$(SIGNATURES)

Check that a covariance matrix is valid for a given `MvNormalLH` object.
"""
function check_cov_matrix(m :: MvNormalLH{T1}, covM :: AbstractMatrix{T1}) where T1
    isValid = true;
    isValid = isValid  && (size(covM) == (n_vars(m), n_vars(m)));
    std2V = sqrt.(diag(covM));
    isValid = isValid  && isapprox(std2V, stds(m), atol = 0.001);
    return isValid
end


"""
	$(SIGNATURES)

Draw random variables given the weight matrix. Returns Matrix by [observation, variable] (unlike Distributions.jl).
"""
function draw_from_weights(m :: MvNormalLH{T1}, 
    wtM :: AbstractMatrix{T1}, nObs :: Integer, rng :: AbstractRNG) where T1

    mvn = make_mv_normal(m, wtM);
    return rand(rng, mvn, nObs)';
end


## -----------  Conditional distribution

"""
    $(SIGNATURES)

Compute conditional distribution of a set of variables, given the others
For Multiple sets of conditioning observations
   
# Arguments
- covM
    covariance matrix
- idx2V
    indices of variables on which we condition
- value2M[observation, variable]
    their values

# Outputs
- condMeanM[observation, variable]
- condStdV(variable)
    conditional means and std of each variable, given all others
    for those not in idx2V
- condCovM(variable, variable)
    conditional covariance; for those not in idx2
"""
function conditional_distrib(m :: MvNormalLH{T1}, 
    idx2V :: AbstractVector{I1}, value2M :: AbstractMatrix{T1}, 
    covM) where {I1 <: Integer, T1}
      
    # Input check
    nObs = size(value2M, 1);
    nCond = length(idx2V);

    @assert all(1 .<= idx2V .<= n_vars(m))
    @assert size(value2M) == (nObs, nCond)
    @assert check_cov_matrix(m, covM)
      
    # Indices of unknowns
    idx1V = [j  for j = 1 : n_vars(m)  if !(j ∈ idx2V)];
      
    sigma11 = covM[idx1V, idx1V];
    sigma12 = covM[idx1V, idx2V];
    sigma21 = covM[idx2V, idx1V];
    sigma22 = covM[idx2V, idx2V];
    mu1V = means(m)[idx1V];
    mu2V = means(m)[idx2V];
    
    condCovM  = sigma11 .- (sigma12 / sigma22) * sigma21;
    condStdV  = sqrt.(diag(condCovM));

    # This is equivalent to mu1 + sigma12 / sigma22 * (value2 - mu2)
    condMeanM = mu1V' .+ (value2M .- mu2V') / sigma22 * sigma21;
      
    # Testing against the original formula; tests matrix expansion
    # This is directly tested by simulation in the test function
    # for i1 = 1 : nObs
    #    condMeanV = (mu1V' .+ ((value2M[i1,:] .- mu2V)' / sigma22) * sigma21)';
    #    @assert (all(abs.(condMeanV .- condMeanM[i1,:]) .< 0.0001));   
    # end
      
    # Output check
    n1 = length(idx1V);
    @assert size(condMeanM) == (nObs, n1)
    @assert size(condStdV) == (n1,)
    @assert size(condCovM) == (n1, n1)

    return condMeanM, condStdV, condCovM
end
 

"""
    $(SIGNATURES)

Weights of conditioning variables for conditional means: `sigma12/sigma22`.
Such that 

    `conditional mean = mu1 + weights * (value2 - mu2)`

# Outputs
- wtM
    n1 x n2 matrix of weights
    each row gives the weights for a "dependent" variable
"""
function cond_mean_weights(m :: MvNormalLH{T1}, 
    idx2V, covM :: AbstractMatrix{T1}) where T1
      
    # Input check
    nVars = size(covM, 1);
    nCond = length(idx2V);

    @assert all(1 .<= idx2V .<= nVars);
    @assert check_cov_matrix(m, covM);

    # Indices of unknowns
    idx1V = [j  for j = 1 : n_vars(m)  if !(j ∈ idx2V)];
      
    sigma12 = covM[idx1V, idx2V];
    sigma22 = covM[idx2V, idx2V];

    # This is equivalent to mu1 + sigma12 / sigma22 * (value2 - mu2)
    wtM = sigma12 / sigma22;
    return wtM
end

# --------------