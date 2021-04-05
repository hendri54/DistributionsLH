Base.@kwdef mutable struct NormalSwitches{T1 <: Real} <: AbstractDistributionSwitches{T1} 
    mean :: T1 = 0.0
    meanLb :: T1 = -10.0
    meanUb :: T1 = 10.0
    meanDescription :: String = "Mean"
    meanLatex :: String = "Mean"
    calMean :: Bool = true
    std :: T1 = 1.0
    stdLb :: T1 = 0.0
    stdUb :: T1 = 10.0
    stdDescription :: String = "Std deviation"
    stdLatex :: String = "Std"
    calStd :: Bool = true
end

mutable struct Normal{T1} <: AbstractDistributionLH{T1}
    objId :: ObjectId
    mean :: T1
    std :: T1
    pvec :: ParamVector
end

## ------------  Generic

draw(u :: Normal{T1}, nDims, rng :: AbstractRNG) where 
    T1 <: AbstractFloat = 
    u.mean .+ u.std .* randn(rng, T1, nDims...);

draw(u :: Normal{T1}, rng :: AbstractRNG)  where T1 <: AbstractFloat = 
    u.mean + u.std * rand(rng, T1);

function quantiles(u :: Normal{T1}, pctV)  where T1 <: AbstractFloat
    normalDistr = Distributions.Normal(u.mean, u.std);
    qV = Distributions.quantile.(normalDistr, pctV);
    return qV
end

Base.eltype(::Normal{T1}) where T1 <: AbstractFloat = T1;

isbounded(u :: Normal{T1}) where T1 = false;

Base.show(io :: IO, um :: Normal{T1}) where T1 = 
    print(io, "Normal{$T1}($(um.mean), $(um.std))");

# Validate draws
function check_draws(u :: Normal{T1}, drawM) where T1
    isValid = true;
    return isValid
end


## ------------  Construction

init_distribution(objId :: ObjectId, switches :: NormalSwitches{T1}) where T1 = 
    init_normal(objId, switches);

function init_normal(objId :: ObjectId, switches :: NormalSwitches{T1}) where T1
    pMean = init_mean(switches);
    pStd = init_std(switches);
    pvec = ParamVector(objId, [pMean, pStd]);
    return Normal(objId, ModelParams.value(pMean), ModelParams.value(pStd), pvec)
end

init_mean(s :: NormalSwitches{T1}) where T1 = 
    Param(:mean, s.meanDescription, s.meanLatex, 
        s.mean, s.mean, s.meanLb, s.meanUb, s.calMean);

init_std(s :: NormalSwitches{T1}) where T1 = 
    Param(:std, s.stdDescription, s.stdLatex, 
        s.std, s.std, s.stdLb, s.stdUb, s.calStd);

make_test_normal_switches() = NormalSwitches(mean = 1.2, std = 2.3);
make_test_normal() = init_normal(ObjectId(:none), make_test_normal_switches());
    
# ---------------