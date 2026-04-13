
"""
    mcsgn(A, method)

Compute the matrix sign of a square matrix `A` using the algorithm `method`.

For a square matrix `A`, the matrix sign `mcsgn(A)` is defined as
```math
mcsgn(A) = (A^2)^{-1/2} A = A (A^2)^{-1/2}
```

If `A` admits an eigendecomposition `A = P Λ P^{-1}`, then `mcsgn(A) = P sign.(Λ) P^{-1}`, i.e. `mcsgn(A)` has the same eigenvectors, but every eigenvalue is replaced by its sign.
"""
function mcsgn end

mcsgn(A::AbstractMatrix, ns::AbstractNewtonSchulz, nsteps=nsteps(ns)) = newton_schulz_square(A, coeffs(ns), nsteps)
