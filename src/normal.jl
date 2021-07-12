mutable struct NormalSwitches{T1} <: AbstractDistributionSwitches{T1} 
    pvec :: ParamVector
    # logNormal :: Bool
end

mutable struct Normal{T1} <: AbstractDistributionLH{T1}
    objId :: ObjectId
    switches :: NormalSwitches{T1}
    mean :: T1
    std :: T1
end


## ------------  Generic

Lazy.@forward Normal.switches (
    ModelParams.get_pvector
);


# function exp_draws!(drawV)
#     for (j, x) in enumerate(drawV)
#         drawV[j] = exp(x);
#     end
# end

function draw(u :: Normal{T1}, nDims, rng :: AbstractRNG) where T1
    drawV = u.mean .+ u.std .* randn(rng, T1, nDims...);
    # is_lognormal(u)  &&  exp_draws!(drawV);
    return drawV
end

function draw(u :: Normal{T1}, rng :: AbstractRNG)  where T1
    drawV = u.mean + u.std * rand(rng, T1);
    # is_lognormal(u)  &&  (drawV = exp(drawV));
    return drawV
end

function quantiles(u :: Normal{T1}, pctV)  where T1
    normalDistr = Distributions.Normal(u.mean, u.std);
    qV = Distributions.quantile.(normalDistr, pctV);
    # is_lognormal(u)  &&  exp_draws!(qV);
    return qV
end

# is_lognormal(switches :: NormalSwitches{T1}) where T1 = switches.logNormal;

Base.eltype(::Normal{T1}) where T1 = T1;

isbounded(u :: Normal{T1}) where T1 = false;

Base.show(io :: IO, um :: Normal{T1}) where T1 = 
    print(io, "Normal{$T1}($(um.mean), $(um.std))");

# Validate draws
function check_draws(u :: Normal{T1}, drawM) where T1
    isValid = true;
    # if is_lognormal(u)
    #     isValid = isValid  &&  all(x -> x > zero(T1), drawM);
    # end
    return isValid
end


## ------------  Construction

ModelParams.get_pvector(switches :: NormalSwitches{T1}) where T1 = switches.pvec;
# ModelParams.get_pvector(nd :: Normal{T1}) where T1 = nd.switches.pvec;

init_distribution(objId :: ObjectId, switches :: NormalSwitches{T1}) where T1 = 
    init_normal(switches);

function init_normal(switches :: NormalSwitches{T1}) where T1
    pvec = get_pvector(switches);
    return Normal{T1}(switches.pvec.objId, switches, 
        param_default_value(pvec, :mean), 
        param_default_value(pvec, :std))
end

function init_normal_switches(objId :: ObjectId; 
    meanVal = 0.0, stdVal = 1.0, calMean = true, calStd = true)
    T1 = typeof(meanVal);
    pMean = init_mean(; meanVal, calMean);
    pStd = init_std(; stdVal, calStd);
    pvec = ParamVector(objId, [pMean, pStd]);
    return NormalSwitches{T1}(pvec)
end

function init_mean(; descr = "Mean", lsym = "μ", meanVal = 0.0,
    meanLb = -10.0, meanUb = 10.0, calMean = true)
    return Param(:mean, descr, lsym, meanVal, meanVal, meanLb, meanUb, calMean);
end

function init_std(; descr = "StdDev", lsym = "σ", stdVal = 1.0,
    stdLb = 0.0, stdUb = 10.0, calStd = true)
    return Param(:std, descr, lsym, stdVal, stdVal, stdLb, stdUb, calStd);
end
        
make_test_normal_switches() = init_normal_switches(ObjectId(:none);
    meanVal = 1.2, stdVal = 2.3);
make_test_normal() = init_normal(make_test_normal_switches());

    
# ---------------