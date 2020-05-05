abstract type AbstractUniform{T1 <: AbstractFloat} <: AbstractDistributionLH end
abstract type AbstractUniformSwitches{T1 <: AbstractFloat} <: AbstractDistributionSwitches end


## -------------------  Generic

draw(u :: AbstractUniform{T1}, nDims, rng :: AbstractRNG) where 
    T1 <: AbstractFloat = 
    min_value(u) .+ value_range(u) .* rand(rng, T1, nDims...);

draw(u :: AbstractUniform{T1}, rng :: AbstractRNG)  where T1 <: AbstractFloat = 
    min_value(u) + value_range(u) * rand(rng, T1);

quantiles(u :: AbstractUniform{T1}, pctV)  where T1 <: AbstractFloat = 
    min_value(u) .+ value_range(u) .* T1.(pctV);

Base.eltype(::AbstractUniform{T1}) where T1 <: AbstractFloat = T1;

Base.show(io :: IO, um :: AbstractUniform{T1}) where T1 = 
    print(io, "Uniform{$T1} over range $(min_value(um)) to $(max_value(um))");

# Validate `Uniform` draws
function check_draws(u :: AbstractUniform{T1}, drawM) where T1
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


## ----------------  Uniform with minimum and range

# The defaults really don't make all that much sense, but we want a keyword constructor
@with_kw mutable struct UniformSwitches{T1 <: AbstractFloat} <: AbstractUniformSwitches{T1}
    xMin :: T1 = 1.0
    xMinLb :: T1 = 0.1
    xMinUb :: T1 = 100.0
    xRange :: T1 = 3.0
    xRangeLb :: T1 = 0.1
    xRangeUb :: T1 = 100.0
    calXMin :: Bool = true
    calXRange :: Bool = true
end


"""
	$(SIGNATURES)

Unform distribution on interval `xMin + [0, xRange]`.
"""
mutable struct Uniform{T1} <: AbstractUniform{T1}
    objId :: ObjectId
    xMin :: T1
    # Difference between h0Max and h0Min
    xRange :: T1
    pvec :: ParamVector
end

min_value(u :: Uniform{T1}) where T1 = u.xMin;
max_value(u :: Uniform{T1}) where T1 = u.xMin + u.xRange;
value_range(u :: Uniform{T1}) where T1 = u.xRange;


## ----------  Construction

function init_uniform(objId :: ObjectId, switches :: UniformSwitches{T1}) where T1
    pMin = init_xmin(switches);
    pRange = init_xrange(switches);
    pvec = ParamVector(objId, [pMin, pRange]);
    xMin = ModelParams.value(pMin);
    return Uniform(objId, xMin, ModelParams.value(pRange), pvec)
end

init_xmin(s :: UniformSwitches{T1}) where T1 = 
    Param(:xMin, "Lower bound", "x_{min}", s.xMin, s.xMin, 
        s.xMinLb, s.xMinUb, s.calXMin);

init_xrange(s :: AbstractUniformSwitches{T1}) where T1 =
    Param(:xRange, "x range", "dx", s.xRange, s.xRange, 
        s.xRangeLb, s.xRangeUb, s.calXRange);

make_test_uniform_switches() = UniformSwitches();
make_test_uniform() = init_uniform(ObjectId(:none),  make_test_uniform_switches());



# With fixed bounds (e.g. HS GPA)
function UniformFixedBounds(objId :: ObjectId, lb :: T1, ub :: T1) where
    T1 <: AbstractFloat

    dx = ub - lb;
    @assert dx > 0
    s = UniformSwitches(
        xMin = lb, xMinLb = lb - 0.1, xMinUb = lb + 0.1, calXMin = false,
        xRange = dx, xRangeLb = dx - 0.1, xRangeUb = dx + 0.1, calXRange = false);
    return init_uniform(objId, s);
    # pMin = Param(:xMin, "Lower bound", "x_{min}", lb, lb, lb - 0.1, lb + 0.1, false);
    # pRange = Param(:xRange, "x range", "dx", dx, dx, dx - 0.1, dx + 0.1, false);
    # pvec = ParamVector(objId, [pMin, pRange])
    # return UniformMarginal(objId, lb, dx, pvec)
end



## ----------------  Uniform with mean and range

@with_kw mutable struct UniformCenteredSwitches{T1 <: AbstractFloat} <: AbstractUniformSwitches{T1}
    xMean :: T1 = 1.0
    xMeanLb :: T1 = 0.1
    xMeanUb :: T1 = 100.0
    xRange :: T1 = 3.0
    xRangeLb :: T1 = 0.1
    xRangeUb :: T1 = 100.0
    calXMean :: Bool = true
    calXRange :: Bool = true
end


"""
	$(SIGNATURES)

Unform distribution on interval `xMean +/- 0.5 * xRange`.
"""
mutable struct UniformCentered{T1} <: AbstractUniform{T1}
    objId :: ObjectId
    xMean :: T1
    # Difference between h0Max and h0Min
    xRange :: T1
    pvec :: ParamVector
end

min_value(u :: UniformCentered{T1}) where T1 = u.xMean - 0.5 * u.xRange;
max_value(u :: UniformCentered{T1}) where T1 = (u.xMean + 0.5 * u.xRange);
value_range(u :: UniformCentered{T1}) where T1 = u.xRange;


## ----------  Construction

function init_uniform(objId :: ObjectId, switches :: UniformCenteredSwitches)
    pMean = init_xmean(switches);
    pRange = init_xrange(switches);
    pvec = ParamVector(objId, [pMean, pRange]);
    xMean = ModelParams.value(pMean);
    return Uniform(objId, xMean, ModelParams.value(pRange), pvec)
end

init_xmean(s :: UniformCenteredSwitches) = 
    Param(:xMean, "Mean", "x_{mean}", s.xMean, s.xMean, 
        s.xMeanLb, s.xMeanUb, s.calXMean);

make_test_uniform_centered_switches() = UniformCenteredSwitches();
make_test_uniform_centered() = init_uniform(ObjectId(:none),  make_test_uniform_centered_switches());



# Count no of observations in quantile groups
# Include 0 as lower bound of first percentile group.
# function count_quantiles(u :: Uniform, drawV, pctV)
#     @assert pctV[1] == 0.0  "Include 0 lower bound"
#     h = fit(Histogram, vec(drawV), quantiles(u, pctV));
#     return h.weights
# end



# """
# 	$(SIGNATURES)

# Converts a N(0,1) input into a Uniform.   
# """
# function scale_normal_endowment(hmS :: Uniform,  inM :: T;
#     checkMarginal :: Bool = true)  where T
    
#     cdfM = cdf.(Normal(), inM);
#     drawM = quantile.(Uniform(hmS.xMin, hmS.xMin + hmS.xRange),  cdfM) :: T
#     return drawM
# end

# -------------