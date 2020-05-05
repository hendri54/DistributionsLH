using Random, Test, ModelParams, DistributionsLH

function uniform_test(um :: AbstractUniform{T1}) where T1 <: AbstractFloat
    println(um)
    T2 = eltype(um);
    @test T1 == T2
    @test isa(max_value(um), T1)

    rng = MersenneTwister(123);

    sizeV = (4, 3);
    uDrawM = draw(um, sizeV, rng);
    @test check_draws(um, uDrawM)
    @test size(uDrawM) == sizeV
    @test eltype(uDrawM) == T1

    uDraw = draw(um, rng);
    @test isa(uDraw, T1)

    qV = quantiles(um, 0.0 : 0.2 : 1.0);
    @test isa(qV, Vector{T1})
    @test isapprox(qV[1], min_value(um))
    @test isapprox(qV[end], max_value(um))
    @test all(diff(qV) .> 0.0)
end


# function um_test()
#     switches = DistributionsLH.make_test_uniform_switches();
#     um = init_uniform(ObjectId(:none), switches);
#     uniform_test(um);
# end

# function um_fixed_bounds_test()
#     um = UniformFixedBounds(ObjectId(:none), 2.0, 4.0);
#     uniform_test(um);
# end

# function centered_test()
#     um = DistributionsLH.make_test_uniform_centered();
#     uniform_test(um);
# end

@testset "Uniform marginal" begin
    for um in [DistributionsLH.make_test_uniform(), 
        DistributionsLH.make_test_uniform_centered(),
        DistributionsLH.make_test_uniform_fixed()]

        uniform_test(um);
    end
end

# ----------