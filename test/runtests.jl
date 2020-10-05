using DistributionsLH, ModelObjectsLH
using Random, StatsBase, GLM, Statistics
using Test

@testset "DistributionsLH" begin
    include("uniform_test.jl");
    include("mv_normal_test.jl");
end

# ------------