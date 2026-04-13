using LinearAlgebra: norm, transpose

"""
    newton_schulz(A, coeffs, nsteps=length(coeffs))

Perform generalized Newton-Schulz iteration on matrix `A`, using the coefficients in `coeffs`.

This method iterates a sequence of quintic polynomials of the form
```math
X_{t+1} = a X_t + b (X_t X_t^T) X_t + c (X_t X_t^T)^2 X_t
```
starting from `X_0 = A / ||A||_F`. 
"""
function newton_schulz(A::AbstractMatrix, coeffs, nsteps = length(coeffs))
    @assert nsteps >= 0
    X = size(A, 1) > size(A, 2) ? transpose(copy(A)) : copy(A)
    X ./= norm(X) + 1e-20
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        Y = X * X'
        X .= a * X + (b * Y + c * (Y*Y)) * X
    end
    return size(A, 1) > size(A, 2) ? transpose(X) : X
end

"""
    newton_schulz!(A, coeffs, nsteps=length(coeffs))

In-place version of `newton_schulz`, which modifies the input matrix `A`.
"""
function newton_schulz!(A::AbstractMatrix, coeffs, nsteps = length(coeffs))
    @assert nsteps >= 0
    X = size(A, 1) > size(A, 2) ? transpose(A) : A
    X ./= norm(X) + 1e-20
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        Y = X * X'
        X .= a * X + (b * Y + c * (Y*Y)) * X
    end
    return size(A, 1) > size(A, 2) ? transpose(X) : X
end

"""
    newton_schulz_square(A, coeffs, nsteps=length(coeffs))

Perform generalized Newton-Schulz iteration on a square matrix `A`, using the coefficients in `coeffs`.

This method iterates a sequence of quintic polynomials of the form
```math
X_{t+1} = a X_t + b  X_t^3 + c X_t^5
```
starting from `X_0 = A / ||A||_F`. 
"""
function newton_schulz_square(A::AbstractMatrix, coeffs, nsteps = length(coeffs))
    @assert nsteps >= 0
    @assert size(A,1) == size(A,2)
    X = copy(A)
    X ./= norm(X) + 1e-20
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        Y = X * X
        X .= a * X + (b * Y + c * (Y*Y)) * X
    end
    return X
end

"""
    newton_schulz_square!(A, coeffs, nsteps=length(coeffs))

In-place version of `newton_schulz_square`, which modifies the input matrix `A`.
"""
function newton_schulz_square!(A::AbstractMatrix, coeffs, nsteps = length(coeffs))
    @assert nsteps >= 0
    @assert size(A,1) == size(A,2)
    X = A
    X ./= norm(X) + 1e-20
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        Y = X * X
        X .= a * X + (b * Y + c * (Y*Y)) * X
    end
    return X
end

function adjust_coeffs(coeffs, factor=1.01)
    return [(a/factor, b/factor^3, c/factor^5) for (a,b,c) in coeffs]
end

abstract type AbstractNewtonSchulz end

nsteps(ns::AbstractNewtonSchulz) = length(coeffs(ns))

"""
    ClassicalNewtonSchulz <: AbstractNewtonSchulz

The classical Newton-Schulz coefficients for approximating the sign function. Convergence is garanteed for any starting value with the unit interval, but is pretty slow.
"""
Base.@kwdef struct ClassicalNewtonSchulz <: AbstractNewtonSchulz
    nsteps::Int = 10
end

const NS_COEFFS = [(2.0, -1.5, 0.5)]

coeffs(::ClassicalNewtonSchulz) = NS_COEFFS
nsteps(ns::ClassicalNewtonSchulz) = ns.nsteps

"""
    NSJordan <: AbstractNewtonSchulz

The coefficients from Keller Jordan's Muon implementation [1]. Convergence is not garranteed, as these coefficients are tuned to bring singular values within the range [0.7,1.3] as fast as possible.

[1] Keller Jordan, *Muon: An optimizer for hidden layers in neural networks*, (2024). https://kellerjordan.github.io/posts/muon/
"""
Base.@kwdef struct NSJordan <: AbstractNewtonSchulz
    nsteps::Int = 5
end

const NSJORDAN_COEFFS = [((3.4445, -4.7750, 2.0315))]

coeffs(::NSJordan) = NSJORDAN_COEFFS
nsteps(ns::NSJordan) = ns.nsteps

"""
    PolarExpress <: AbstractNewtonSchulz

The Polar Express coefficients from [1]. These coefficients are designed to approximate the matrix sign function with precise control of the absolute error on singular values.

[1] N. Amsel, D. Persson, C. Musco and R. M. Gower, *The Polar Express: Optimal Matrix Sign Methods and Their Application to the Muon Algorithm*, (2025). https://arxiv.org/abs/2505.16932
"""
struct PolarExpress <: AbstractNewtonSchulz end

const POLAREXPRESS_COEFFS = adjust_coeffs([
    (8.28721201814563 , -23.595886519098837 , 17.300387312530933) ,
    (4.107059111542203 , -2.9478499167379106 , 0.5448431082926601) ,
    (3.9486908534822946 , -2.908902115962949 , 0.5518191394370137) ,
    (3.3184196573706015 , -2.488488024314874 , 0.51004894012372) ,
    (2.300652019954817 , -1.6689039845747493 , 0.4188073119525673) ,
    (1.891301407787398 , -1.2679958271945868 , 0.37680408948524835) ,
    (1.8750014808534479 , -1.2500016453999487 , 0.3750001645474248) ,
    (1.875 * 1.01, -1.25 * 1.01^2, 0.375 * 1.01^3) # don't adjust the last coefficients
])

coeffs(::PolarExpress) = POLAREXPRESS_COEFFS

"""
    NSJianlinSu <: AbstractNewtonSchulz

The Newton-Schulz coefficients from Jianlin Su [1]. These coefficients are obtained by the same method as the PolarExpress coefficients, and are indeed very close.

[1] Jianlin Su, *Newton-Schulz Iteration for the msign Operator (Part 2), (2025). https://kexue.fm/archives/10996, english translation: https://datasets.osmarks.net/kexue/site/10996-Newton-Schulz-Iteration-for-the-msign-Operator-Part-2.html
"""
struct NSJianlinSu <: AbstractNewtonSchulz end

const NSJIANLINSU_COEFFS = adjust_coeffs([
    (8.28721, -23.5959, 17.3004),
    (4.10706, -2.94785, 0.544843),
    (3.94869, -2.9089, 0.551819),
    (3.31842, -2.48849, 0.510049),
    (2.30065, -1.6689, 0.418807),
    (1.8913, -1.268, 0.376804),
    (1.875, -1.25, 0.375),
    (1.875, -1.25, 0.375)
])

coeffs(::NSJianlinSu) = NSJIANLINSU_COEFFS


struct NSCesista <: AbstractNewtonSchulz end

const NSCESISTA_COEFFS = [
    (7.2086, -15.5131, 9.0178),
    (3.9623, -2.5813, 0.4542),
    (3.9466, -2.5765, 0.4544),
    (3.8991, -2.5671, 0.4566),
    (3.7186, -2.5308, 0.4653),
    (3.1390, -2.3073, 0.4733),
    (2.1715, -1.5246, 0.3885),
    (1.8648, -1.2224, 0.3577),
]

coeffs(::NSCesista) = NSCESISTA_COEFFS

struct NSSquareRoot <: AbstractNewtonSchulz end

const NSSQRT_COEFFS = [
    (7.424865680309214, -18.39581635618996, 12.896720413604342),
    (3.4877256051546017, -2.3300436563986993, 0.4404692168431095),
    (2.7766085124882527, -2.070643152532662, 0.46302261050004967),
    (1.9913142104341506, -1.373936700681269, 0.3875934979568538),
    (1.8754637749479246, -1.2505152090010534, 0.37505152463617264),
    (1.874999066623701, -1.2499981332141676, 0.37499906659046633),
    (1.875, -1.25, 0.375),
]

coeffs(::NSSquareRoot) = NSSQRT_COEFFS
