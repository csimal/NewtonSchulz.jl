using LinearAlgebra: eigen, norm, opnorm, Diagonal
@testset "mcsgn" begin

    @testset "msign Newton-Schulz" begin
        rng = Xoshiro(2026)
        A = randn(10, 10)
        A = A + A' # make A symmetric
        F = eigen(A)
        B = F.vectors * Diagonal(sign.(F.values)) * F.vectors'

        @testset "Classical Newton-Schulz" begin
            C = mcsgn(A, ClassicalNewtonSchulz())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 5e-03
            @test opnorm(B-C) < 1e-02
        end
        @testset "Jordan Newton-Schulz" begin
            C = mcsgn(A, NSJordan())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 0.5
            @test opnorm(B-C) < 0.5
        end
        @testset "PolarExpress" begin
            C = mcsgn(A, PolarExpress())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 5e-03 # surprisingly poor absolute tolerance
            @test opnorm(B-C) < 5e-03
        end
        @testset "JianlinSu" begin
            C = mcsgn(A, NSJianlinSu())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 1e-03
            @test opnorm(B-C) < 1e-03
        end
        @testset "Cesista" begin
            C = mcsgn(A, NSCesista())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 1e-03
            @test opnorm(B-C) < 1e-03
        end
        @testset "Square Root" begin
            C = mcsgn(A, NSSquareRoot())

            @test size(C) == size(A)
            @test norm(B-C, Inf) < 5e-02
            @test opnorm(B-C) < 5e-02
        end
    end
end
