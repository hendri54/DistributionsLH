mutable struct DegenerateSwitches{T1 <: AbstractFloat} <: AbstractDistributionSwitches{T1}
    fixedValue :: T1
end

mutable struct DegenerateDistribution{T1 <: AbstractFloat} <: AbstractDistributionLH{T1}
    objId :: ObjectId
    fixedValue :: T1
end

ModelParams.has_pvector(u :: DegenerateDistribution) = false;

draw(u :: DegenerateDistribution, nDims, rng :: AbstractRNG) = 
    fill(u.fixedValue, nDims...);

draw(u :: DegenerateDistribution, rng :: AbstractRNG) = u.fixedValue;

max_value(u :: DegenerateDistribution) = u.fixedValue;
min_value(u :: DegenerateDistribution) = u.fixedValue;

quantiles(u :: DegenerateDistribution{T1}, pctV)  where T1 <: AbstractFloat = 
    fill(u.fixedValue, size(pctV));

Base.eltype(::DegenerateDistribution{T1}) where T1 <: AbstractFloat = T1;

isbounded(u :: DegenerateDistribution{T1}) where T1 = true;

Base.show(io :: IO, um :: DegenerateDistribution{T1}) where T1 = 
    print(io, "Degenerate{$T1} with fixed value $(um.fixedValue)");

function describe(switches :: DegenerateSwitches)
    return ["Degenerate distribution"  "fixed value $(switches.fixedValue)"];
end

function describe(switches :: DegenerateDistribution)
    return ["Degenerate distribution"  "fixed value $(switches.fixedValue)"];
end

# Validate `Uniform` draws
function check_draws(u :: DegenerateDistribution{T1}, drawM) where T1
    isValid = true;
    for draw in drawM
        isapprox(draw, u.fixedValue)  ||  (isValid = false);
    end
    return isValid
end

function init_distribution(objId :: ObjectId, switches :: DegenerateSwitches) 
    return DegenerateDistribution(objId, switches.fixedValue);
end

make_test_degenerate_switches() = DegenerateSwitches(1.0);
make_test_degenerate() = init_distribution(make_test_degenerate_switches());


# ---------------