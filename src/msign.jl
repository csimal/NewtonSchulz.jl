
"""
    msign(A, method)

Compute the matrix sign of a rectangular matrix `A` using the algorithm `method`.

For a rectangular matrix `A`, `msign(A)` is defined as
```math
msign(A) = (A * A^T)^{-1/2} A = A (A * A^T)^{-1/2} 
```
If `A = U * Σ * V^T` is the SVD decomposition of `A`, then `msign(A) = U sign.(Σ) * V^T`, i.e. `msign(A)` has the same singular vectors as `A` but any non-zero singular value is replaced by one.
"""
function msign end

struct MSignSVD end

function msign(A::AbstractMatrix, ::MSignSVD)
    F = svd(A) # computes thin SVD by default
    return F.U * F.Vt
end

msign(A::AbstractMatrix, ns::AbstractNewtonSchulz, nsteps=nsteps(ns)) = newton_schulz(A, coeffs(ns), nsteps)
