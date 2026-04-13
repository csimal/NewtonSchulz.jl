@testset "Newton-Schulz" begin
    @testset "adjust_coeffs" begin
        using NewtonSchulz: adjust_coeffs
        coeffs = [(1.0,1.0,1.0), (2.0,2.0,2.0)]
        adjusted_coeffs = adjust_coeffs(coeffs)
        @test length(adjusted_coeffs) == length(coeffs)
        @test eltype(adjusted_coeffs) == eltype(coeffs)
        @test adjusted_coeffs[1] == (1/1.01, 1/1.01^3, 1/1.01^5)
        @test adjusted_coeffs[2] == (2/1.01, 2/1.01^3, 2/1.01^5)

        adjusted_coeffs = adjust_coeffs(coeffs, 1.0)
        @test all(adjusted_coeffs .== coeffs)
    end
    @testset "newton_schulz" begin
        using NewtonSchulz: NS_COEFFS
        coeffs = NS_COEFFS
        rng = Xoshiro(2026)

        A = rand(rng, 10, 20)
        B = newton_schulz(A, coeffs)

        @test size(B) == size(A)
        @test B isa Matrix

        B = newton_schulz(A', coeffs)
        @test (size(B,2), size(B,1)) == size(A)

        B = newton_schulz(A, coeffs, 10)
        @test size(B) == size(A)

        A = rand(20, 10)
        B = newton_schulz(A, coeffs)
        @test size(B) == size(A)

        @test_throws AssertionError newton_schulz(A, coeffs, -1)
    end
    @testset "newton_schulz!" begin
        using NewtonSchulz: NS_COEFFS
        coeffs = NS_COEFFS
        rng = Xoshiro(2026)

        A = rand(rng, 10, 20)
        B = newton_schulz!(A, coeffs)

        @test size(B) == size(A)
        @test B isa Matrix
        @test A === B

        B = newton_schulz!(A', coeffs)
        @test (size(B,2), size(B,1)) == size(A)

        C = newton_schulz!(A, coeffs, 10)
        @test size(C) == size(A)

        A = rand(20, 10)
        B = newton_schulz!(A, coeffs)
        @test size(B) == size(A)
        @test A === B

        @test_throws AssertionError newton_schulz!(A, coeffs, -1)
    end
    @testset "newton_schulz_square" begin
        using NewtonSchulz: NS_COEFFS, newton_schulz_square
        coeffs = NS_COEFFS
        rng = Xoshiro(2026)
        A = randn(rng, 10,10)
        B = newton_schulz_square(A, coeffs)

        @test size(B) == size(A)

        B = newton_schulz_square(A, coeffs, 10)

        @test size(B) == size(A)

        @test_throws AssertionError newton_schulz_square(A, coeffs, -1)
    end
    @testset "newton_schulz_square!" begin
        using NewtonSchulz: NS_COEFFS, newton_schulz_square
        coeffs = NS_COEFFS
        rng = Xoshiro(2026)
        A = randn(rng, 10,10)
        B = newton_schulz_square!(A, coeffs)

        @test size(B) == size(A)
        @test B === A

        A = randn(rng, 10, 10)
        C = newton_schulz_square!(A, coeffs, 10)

        @test size(C) == size(A)
        @test C === A

        @test_throws AssertionError newton_schulz_square!(A, coeffs, -1)
    end
end
