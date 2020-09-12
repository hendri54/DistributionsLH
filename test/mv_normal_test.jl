using DistributionsLH
using GLM, StatsBase, Random, Test

d = DistributionsLH;

function pretty_print(x)
    show(stdout, MIME("text/plain"), x);
    println("");
end

# The numbers here need to be nice. Otherwise it takes very large samples to get accurate cov matrices.
function make_mean_std(n :: Integer)
    meanV = collect(range(-2.0, 3.0, length = n));
    stdV  = collect(range(0.4, 1.2, length = n));
    # if length(stdV) > 2
    #     stdV[2] = 0.0;
    # end
    return meanV, stdV
end


# Make lower triangular weight matrix
function make_weight_matrix(stdV)
   n = length(stdV);
   rng = MersenneTwister(123);
   wtM = randn(rng, n,n) .- 0.5;
   for i1 = 1 : n
      if i1 < n
         wtM[i1, (i1+1) : n] .= 0.0;
      end
      if stdV[i1] == 0.0
         wtM[i1, :] .= 0.0;
         wtM[:, i1] .= 0.0;
      end
      # Diagonal elements are 1
      wtM[i1, i1] = 1.0;
   end
   return wtM
end


function mvn_basic_test(n :: Integer)
    @testset "Basics" begin
        meanV, stdV = make_mean_std(n);
        switches = MvNormalLHSwitches(meanV, stdV);
        m = MvNormalLH(switches);
        @test isapprox(d.means(m), meanV)
        println(m);

        wtM = make_weight_matrix(stdV);
        @test check_weight_matrix(m, wtM)
        covM = cov_matrix(m, wtM);
        @test check_cov_matrix(m, covM)
    end
end


function drawing_test(n :: Integer)
    meanV, stdV = make_mean_std(n);
    switches = MvNormalLHSwitches(meanV, stdV);
    m = MvNormalLH(switches);
    rng = MersenneTwister(43);
    nObs = 50_000;
    wtM = make_weight_matrix(stdV);
    drawM = draw_from_weights(m, wtM, nObs, rng);
    @test isequal(size(drawM), (nObs, n_vars(m)))

    mean2V = vec(mean(drawM, dims = 1));
    @test isapprox(mean2V, meanV, atol = 0.03)
    std2V = vec(StatsBase.std(drawM, dims = 1))
    @test isapprox(std2V, stdV, rtol = 0.02)

    covM = Statistics.cov(drawM, dims = 1);
    cov2M = cov_matrix(m, wtM);

    maxDiff = maximum((abs.(covM .- cov2M) ./ max.(0.2, cov2M)));
    @test maxDiff < 0.05
    if maxDiff > 0.05
        pretty_print(covM)
        pretty_print(cov2M)
        pretty_print((covM .- cov2M) ./ max.(0.2, cov2M))
    end
end


function conditional_test()
    @testset "Conditional distribution" begin
        n = 5;
        meanV, stdV = make_mean_std(n);
        switches = MvNormalLHSwitches(meanV, stdV);
        m = MvNormalLH(switches);
        wtM = make_weight_matrix(stdV);
        covM = cov_matrix(m, wtM);

        idx2V = [2, 4];
        value2V = [0.5 -0.3];
        idx1V = [j  for j = 1 : n  if !(j ∈ idx2V)];
        n1 = length(idx1V);
        conditional_distrib(m, idx2V, value2V, covM);
     
        # Draw a sample, by [variable, observation]
        rng = MersenneTwister(43);
        nObs = 10_000;
        # Transpose so we have draws by [obs, var]
        drawM = draw_from_weights(m, wtM, nObs, rng);
        
        condMeanM = zeros(nObs, n1);
        condStdM  = zeros(nObs, n1);
        for i1 = 1 : nObs
            # Columns are variables
            value2V = drawM[i1 : i1, idx2V];
            condMeanM[i1, :], condStdM[i1,:], _ = 
                conditional_distrib(m, idx2V, value2V, covM);
        end
        
        # Test with multiple inputs (tests matrix expansion)
        condMean2M, _ = conditional_distrib(m, idx2V, drawM[:, idx2V], covM);
        @test isapprox(condMean2M, condMeanM,  atol = 0.0001);
     
        # ----  Check by simulation + regression
        for iVar = 1 : n1
            # Regress y against its conditional mean
            idx1 = idx1V[iVar];
            linModel = lm(hcat(ones(nObs), condMeanM[:, iVar]),  drawM[:, idx1]);
            betaV = coef(linModel);
            seBetaV = stderror(linModel);
            
            # Test that intercept is 0 and slope is 1
            @test !any(abs.(betaV .- [0.0, 1.0]) .> 2.0 * seBetaV)
            
            # Residuals
            residV = drawM[:, idx1] .- condMeanM[:, iVar];
            std2 = std(residV);
            @test abs(std2 - mean(condStdM[:, iVar])) < 0.01
        end         
	end
end


# Conditional mean weights
function cond_mean_weights_test()
    @testset "Conditional mean weights" begin
        rng = MersenneTwister(42);
        n = 7;
        # With means of 0 and conditioning values of 1, we can compare
        # conditional means and weights for conditional means
        meanV = zeros(n);
        stdV  = rand(n) * 2.0;
        
        m = MvNormalLH(meanV, stdV);
        wtM = make_weight_matrix(stdV);
        covM = cov_matrix(m, wtM);

        # Conditioning indices
        idx2V = [2, 4, 5];
        value2V = ones(1, length(idx2V));
        condMeanV, _ = conditional_distrib(m, idx2V, value2V, covM);
        # Compare with conditioning weights
        wtM = cond_mean_weights(m, idx2V, covM);
        
        @test isapprox(vec(condMeanV), vec(sum(wtM, dims = 2)), atol = 1e-4);
    end
end



@testset "MvNormal" begin
    for n ∈ [2, 5]
        mvn_basic_test(n);
        drawing_test(n);
    end
    conditional_test();
    cond_mean_weights_test();
end

# ---------------