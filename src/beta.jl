abstract type AbstractBeta{T1} <: AbstractDistributionLH{T1} end
abstract type AbstractBetaSwitches{T1} <: AbstractDistributionSwitches{T1} end


## -------------------  Generic

function draw(u :: AbstractBeta{T1}, nDims, rng :: AbstractRNG) where 
    T1 <: AbstractFloat
    betaDistr = Distributions.Beta(u.alpha, u.beta);
    draws = min_value(u) .+ value_range(u) .* rand(rng, betaDistr, nDims...);
    return draws
end

function draw(u :: AbstractBeta{T1}, rng :: AbstractRNG)  where T1 <: AbstractFloat 
    betaDistr = Distributions.Beta(u.alpha, u.beta);
    draws = min_value(u) + value_range(u) * rand(rng, betaDistr);
    return draws
end

function quantiles(u :: AbstractBeta{T1}, pctV)  where T1 <: AbstractFloat
    betaDistr = Distributions.Beta(u.alpha, u.beta);
    qV = min_value(u) .+ value_range(u) .* Distributions.quantile.(betaDistr, pctV);
    return qV
end

isbounded(u :: AbstractBeta{T1}) where T1 = true;
min_value(u :: AbstractBeta{T1}) where T1 = u.xMin;
max_value(u :: AbstractBeta{T1}) where T1 = u.xMin + u.xRange;
value_range(u :: AbstractBeta{T1}) where T1 = u.xRange;
alpha_value(u :: AbstractBeta{T1}) where T1 = u.alpha;
beta_value(u :: AbstractBeta{T1}) where T1 = u.beta;

Base.eltype(::AbstractBeta{T1}) where T1 <: AbstractFloat = T1;

Base.show(io :: IO, um :: AbstractBeta{T1}) where T1 = 
    print(io, "Beta{$T1} over range $(min_value(um)) to $(max_value(um))");

# Validate draws
function check_draws(u :: AbstractBeta{T1}, drawM) where T1
    isValid = true;
    if any(x -> x < min_value(u), drawM)
        @warn "Below $(min_value(u)):  $(minimum(drawM))"
        isValid = false;
    end

    if any(x -> x > max_value(u), drawM)
        @warn "Above $(max_value(u)):  $(maximum(drawM))"
        isValid = false;
    end
    return isValid
end


init_distribution(objId :: ObjectId, switches :: AbstractBetaSwitches{T1}) where T1 =
    init_beta(objId, switches);

make_test_beta_switches() = BetaSwitches();
make_test_beta() = init_beta(ObjectId(:none),  make_test_beta_switches());
    

## ----------  Default Beta

"""
	$(SIGNATURES)

Switches for default Beta distribution object.
"""
Base.@kwdef mutable struct BetaSwitches{T1 <: AbstractFloat} <: AbstractBetaSwitches{T1}
    xMin :: T1 = 1.0
    xMinLb :: T1 = 0.1
    xMinUb :: T1 = 100.0
    xMinDescription :: String = "Lower bound"
    xMinLatex :: String = "x_{min}"
    xRange :: T1 = 3.0
    xRangeLb :: T1 = 0.1
    xRangeUb :: T1 = 100.0
    xRangeDescription :: String = "Range"
    xRangeLatex :: String = "dx"
    alpha :: T1 = 2.0
    alphaLb :: T1 = 1.0
    alphaUb :: T1 = 3.0
    alphaDescription :: String = "alpha"
    alphaLatex :: String = "\\alpha"
    beta :: T1 = 2.0
    betaLb :: T1 = 1.0
    betaUb :: T1 = 5.0
    betaDescription :: String = "beta"
    betaLatex :: String = "\\beta"
    calXMin :: Bool = true
    calXRange :: Bool = true
    calAlpha :: Bool = true
    calBeta :: Bool = true
end


"""
	$(SIGNATURES)

Beta marginal distribution. Characterized by lower bound, upper bound, and Beta parameters alpha and beta.
Cannot store the Beta distribution in the object because its parameters are not known.
"""
mutable struct Beta{T1} <: AbstractBeta{T1}
    objId :: ObjectId
    pvec :: ParamVector
	xMin :: T1
    xRange :: T1
    alpha :: T1
    beta :: T1
	# betaDistr :: Distributions.Beta{T1}
end


"""
	$(SIGNATURES)

Construct a Uniform distribution from its switches.
"""    
function init_beta(objId :: ObjectId, switches :: BetaSwitches{T1}) where T1
    pMin = init_xmin(switches);
    pRange = init_xrange(switches);
    pAlpha = init_alpha_param(switches);
    pBeta = init_beta_param(switches);
    pvec = ParamVector(objId, [pMin, pRange, pAlpha, pBeta]);
    xMin = ModelParams.value(pMin);
    # betaDistr = Distributions.Beta(alpha, beta);
    return Beta(objId, pvec, 
        ModelParams.value(pMin), ModelParams.value(pRange), 
        ModelParams.value(pAlpha), ModelParams.value(pBeta));
end

init_xmin(s :: AbstractBetaSwitches{T1}) where T1 = 
    Param(:xMin, s.xMinDescription, s.xMinLatex, s.xMin, s.xMin, 
        s.xMinLb, s.xMinUb, s.calXMin);

init_xrange(s :: AbstractBetaSwitches{T1}) where T1 =
    Param(:xRange, s.xRangeDescription, s.xRangeLatex, s.xRange, s.xRange, 
        s.xRangeLb, s.xRangeUb, s.calXRange);

init_alpha_param(s :: AbstractBetaSwitches{T1}) where T1 =
    Param(:alpha, s.alphaDescription, s.alphaLatex, s.alpha, s.alpha, 
        s.alphaLb, s.alphaUb, s.calAlpha);

init_beta_param(s :: AbstractBetaSwitches{T1}) where T1 =
    Param(:beta, s.betaDescription, s.betaLatex, s.beta, s.beta, 
        s.betaLb, s.betaUb, s.calBeta);
            


# --------------