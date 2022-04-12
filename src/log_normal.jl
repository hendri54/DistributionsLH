mutable struct LogNormalSwitches{T1} <: AbstractDistributionSwitches{T1} 
    pvec :: ParamVector
end

mutable struct LogNormal{T1} <: AbstractDistributionLH{T1}
    objId :: ObjectId
    switches :: LogNormalSwitches{T1}
    lb :: T1
    mean :: T1
    std :: T1
end

Lazy.@forward LogNormal.switches (
    ModelParams.get_pvector
);


function draw(u :: LogNormal{T1}, nDims, rng :: AbstractRNG) where T1
    drawV = u.lb .+ exp.(u.mean .+ u.std .* randn(rng, T1, nDims...));
    # is_lognormal(u)  &&  exp_draws!(drawV);
    return drawV
end

function draw(u :: LogNormal{T1}, rng :: AbstractRNG)  where T1
    drawV = u.lb + exp(u.mean + u.std * rand(rng, T1));
    # is_lognormal(u)  &&  (drawV = exp(drawV));
    return drawV
end

function quantiles(u :: LogNormal{T1}, pctV)  where T1
    normalDistr = Distributions.Normal(u.mean, u.std);
    qV = u.lb .+ exp.(Distributions.quantile.(normalDistr, pctV));
    # is_lognormal(u)  &&  exp_draws!(qV);
    return qV
end

Base.eltype(::LogNormal{T1}) where T1 = T1;

isbounded(u :: LogNormal{T1}) where T1 = false;

Base.show(io :: IO, um :: LogNormal{T1}) where T1 = 
    print(io, "LogNormal{$T1}($(um.mean), $(um.std))");

# Validate draws
function check_draws(u :: LogNormal{T1}, drawM) where T1
    isValid = true;
    isValid = isValid  &&  all(x -> x >= u.lb, drawM);
    return isValid
end


## ----------  Construction

ModelParams.get_pvector(switches :: LogNormalSwitches{T1}) where T1 = switches.pvec;

init_distribution(objId :: ObjectId, switches :: LogNormalSwitches{T1}) where T1 = 
    init_log_normal(switches);

function init_log_normal(switches :: LogNormalSwitches{T1}) where T1
    pvec = get_pvector(switches);
    return LogNormal{T1}(switches.pvec.objId, switches, 
        param_default_value(pvec, :lb),
        param_default_value(pvec, :mean), 
        param_default_value(pvec, :std))
end

function init_log_normal_switches(objId :: ObjectId; 
    lbVal = 1.0, calLb = false, 
    meanVal = 0.0, calMean = true, 
    stdVal = 0.5, stdLb = 0.1, stdUb = 1.0, calStd = true
    )
    T1 = typeof(meanVal);
    pLb = init_lb(; lbVal, calLb);
    pMean = init_mean(; meanVal, calMean);
    pStd = init_std(; stdVal, stdLb, stdUb, calStd);
    pvec = ParamVector(objId, [pLb, pMean, pStd]);
    return LogNormalSwitches{T1}(pvec)
end

function init_lb(; descr = "Lower bound", lsym = "lb", lbVal = 0.0,
    lbLb = -10.0, lbUb = 10.0, calLb = false)
    return Param(:lb, descr, lsym, lbVal, lbVal, lbLb, lbUb, calLb);
end

make_test_log_normal_switches() = init_log_normal_switches(ObjectId(:none);
    meanVal = 1.2, stdVal = 2.3); 
make_test_log_normal() = init_log_normal(make_test_log_normal_switches());


# ---------------------