# Source: https://datasets.osmarks.net/kexue/site/11158-Efficient-Computation-of-Matrix-Square-Roots-and-Inverse-Square-Roots.html

"""
    msqrt(A, method)

Compute the matrix square root of a positive semidefinite square matrix `A`.

The square root of a matrix `A` is a matrix `B` such that `A = B^2`. This function computes the unique square root whose eigenvalues are all positive (i.e. the arithmetic square root)
"""
function msqrt end

msqrt(A::AbstractMatrix, ns::AbstractNewtonSchulz) = msqrt(A, coeffs(ns))
msqrt(A::AbstractMatrix, ns::ClassicalNewtonSchulz) = msqrt(A, coeffs(ns), ns.nsteps)
msqrt(A::AbstractMatrix, ns::NSJordan) = msqrt(A, coeffs(ns), ns.nsteps)

function msqrt(A::AbstractMatrix, coeffs, nsteps=length(coeffs))
    @assert nsteps >= 0
    @assert size(A,1) == size(A,2) "matrix must be square"
    T = tr(A)
    Y = copy(A) / T
    YZ = copy(Y)
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        W = a * I + b * YZ + c * YZ * YZ
        Y = W * Y
        YZ = W * W * YZ
    end
    return Y * sqrt(T)
end


"""
    mrsqrt(A, method)
    
Compute the inverse square root of a positive semidefinite matrix `A`.
"""
function mrsqrt end

mrsqrt(A::AbstractMatrix, ns::AbstractNewtonSchulz) = msqrt(A, coeffs(ns))
mrsqrt(A::AbstractMatrix, ns::ClassicalNewtonSchulz) = msqrt(A, coeffs(ns), ns.nsteps)
mrsqrt(A::AbstractMatrix, ns::NSJordan) = msqrt(A, coeffs(ns), ns.nsteps)

function mrsqrt(A::AbstractMatrix, coeffs, nsteps=length(coeffs))
    @assert nsteps >= 0
    @assert size(A,1) == size(A,2) "matrix must be square"
    T = tr(A)
    YZ = A / T
    Z = zeros(size(A)) + I
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        W = a * I + b * YZ + c * YZ * YZ
        Z = Z * W
        YZ = W * W * YZ
    end
    return Z / sqrt(T)
end

"""
    matmul_mrsqrt(A, B, method)

Compute the product `A * B^{-1/2}` where `A` is a rectangular matrix, and `B` is a positive definite matrix.
"""
function matmul_mrsqrt end

matmul_mrsqrt(A, B, ns::AbstractNewtonSchulz) = matmul_mrsqrt(A, B, coeffs(ns))
matmul_mrsqrt(A, B, ns::ClassicalNewtonSchulz) = matmul_mrsqrt(A, B, coeffs(ns), ns.nsteps)
matmul_mrsqrt(A, B, ns::NSJordan) = matmul_mrsqrt(A, B, coeffs(ns), ns.nsteps)

function matmul_mrsqrt(A, B, coeffs, nsteps=length(coeffs))
    @assert size(A,2) == size(B,1)
    @assert size(B,1) == size(B,2)
    T = tr(B)
    YZ = B / T
    Z = copy(A)
    for t in 1:nsteps
        (a,b,c) = t > length(coeffs) ? coeffs[end] : coeffs[t]
        W = a * I + b * YZ + c * YZ * YZ
        Z = Z * W
        YZ = W * W * YZ
    end
    return Z / sqrt(T)
end
