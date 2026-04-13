using LinearAlgebra: svd, norm, opnorm
@testset "msign" begin
    @testset "msign SVD" begin
        rng = Xoshiro(2026)
        A = randn(rng, 10, 20)
        B = msign(A, MSignSVD())

        @test size(B) == size(A)
        @test B isa Matrix

        F_A = svd(A)
        F_B = svd(B)
        @test all(isapprox.(F_B.S, 1.0, atol=0.001))
    end
    @testset "msign Newton-Schulz" begin
        rng = Xoshiro(2026)
        A = randn(10, 20)
        B = msign(A, MSignSVD())

        @testset "Classical Newton-Schulz" begin
            C = msign(A, ClassicalNewtonSchulz())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 1e-03
            @test opnorm(B-C) < 1e-03
        end
        @testset "Jordan Newton-Schulz" begin
            C = msign(A, NSJordan())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 0.5
            @test opnorm(B-C) < 0.5
        end
        @testset "PolarExpress" begin
            C = msign(A, PolarExpress())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 5e-03 # surprisingly poor absolute tolerance
            @test opnorm(B-C) < 5e-03
        end
        @testset "JianlinSu" begin
            C = msign(A, NSJianlinSu())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 1e-03
            @test opnorm(B-C) < 1e-03
        end
        @testset "Cesista" begin
            C = msign(A, NSCesista())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 1e-03
            @test opnorm(B-C) < 1e-03
        end
    end
end
