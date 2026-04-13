using LinearAlgebra: norm, opnorm, I
using NewtonSchulz: coeffs

@testset "matrix square root" begin
    rng = Xoshiro(2026)
    A = randn(rng, 10, 10) ./ 10
    A = A * A' + I # ensure symmetric positive definite

    @testset "msqrt" begin
        B = msqrt(A, NSJianlinSu())

        @test norm(A - B^2, Inf) < 1e-03
        @test opnorm(A - B^2) < 1e-3

        coeffs_ = coeffs(NSJianlinSu())
        C = msqrt(A, coeffs_, 10)

        @test norm(A - C^2, Inf) < 1e-03
        @test opnorm(A - C^2) < 1e-03

        @test_throws AssertionError msqrt(A, coeffs_, -1)
        @test_throws AssertionError msqrt(ones(3,4), coeffs_, 10)
    end
    @testset "mrsqrt" begin
        B = mrsqrt(A, NSJianlinSu())

        # terrible tolerances
        @test norm(I - A * B^2, Inf) < 1.0
        @test opnorm(I - A * B^2) < 1.0

        B = mrsqrt(A, NSSquareRoot())
        @test norm(I - A * B^2, Inf) < 1.0
        @test opnorm(I - A * B^2) < 1.0

        coeffs_ = coeffs(NSJianlinSu())
        C = mrsqrt(A, coeffs_, 10)

        @test norm(I - A * C^2, Inf) < 1e-03
        @test opnorm(I - A * C^2) < 1e-03

        @test_throws AssertionError mrsqrt(A, coeffs_, -1)
        @test_throws AssertionError mrsqrt(ones(3,4), coeffs_, 10)
    end
    @testset "matmul_mrsqrt" begin

    end
end
