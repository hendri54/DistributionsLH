using Random, Test, ModelParams, DistributionsLH

mdl = DistributionsLH;

function distrib_test(switches :: AbstractDistributionSwitches{T1}) where T1
    um = init_distribution(ObjectId(:distrib), switches);
    println(um)
    T2 = eltype(um);
    @test T1 == T2

    if isbounded(um)
        @test isa(max_value(um), T1)
        @test min_value(um) < max_value(um)
    end

    rng = MersenneTwister(123);

    sizeV = (4, 3);
    uDrawM = draw(um, sizeV, rng);
    @test check_draws(um, uDrawM)
    @test size(uDrawM) == sizeV
    @test eltype(uDrawM) == T1

    uDraw = draw(um, rng);
    @test isa(uDraw, T1)

    if isbounded(um)
        qV = quantiles(um, 0.0 : 0.2 : 1.0);
        @test isa(qV, Vector{T1})
        @test isapprox(qV[1], min_value(um))
        @test isapprox(qV[end], max_value(um))
        @test all(diff(qV) .> 0.0)
    end

    inM = randn(rng, sizeV...);
    drawM = scale_normal_draws(um, inM);
    @test size(inM) == size(drawM)
    @test check_draws(um, drawM)
end


@testset "Uniform marginal" begin
    for switches in [mdl.make_test_uniform_switches(), 
        mdl.make_test_uniform_centered_switches()]

        distrib_test(switches);
    end
end

@testset "Beta" begin
    for switches in [mdl.make_test_beta_switches()]
        distrib_test(switches);
    end
end

@testset "Normal" begin
    for switches in (mdl.make_test_normal_switches(),
        mdl.make_test_log_normal_switches())
        distrib_test(switches);
    end
end

# ----------