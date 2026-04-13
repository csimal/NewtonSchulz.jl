# NewtonSchulz

[![Build Status](https://github.com/csimal/NewtonSchulz.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/csimal/NewtonSchulz.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/csimal/NewtonSchulz.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/csimal/NewtonSchulz.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This package provides a generic implementation of the Newton-Schulz method for computing various matrix functions such as the matrix sign and matrix square root.

## Feature Summary

This package provides iterative methods for computing various matrix functions, for example, the following computes the matrix sign of a rectangular matrix `A` using the PolarExpress coefficients.

```julia
O = msign(A, PolarExpress())
```

All provided methods perform a default number of iterations, which can be over-riden as follows

```julia
O = msign(A, PolarExpress(), nsteps)
```

All the matrix functions in this package follow this interface.

### Functions exported by this package

- `msign`: Compute the matrix sign of a rectangular matrix
- `mcsgn`: Compute the matrix sign of a square matrix
- `msqrt`: Compute the arithmetic square root of a square positive semidefinite matrix
- `mrsqrt`: Compute the inverse square root of a square positive semidefinite matrix
- `matmul_mrsqrt`: Compute the product $A B^{-1/2}$ of a rectangular matrix $A$ and a square positive semidefinite matrix $B$
- `newton_schulz` and `newton_schulz!`: Perform a generic Newton-Schulz iteration
- `newton_schulz_square` and `newton_schulz_square!`: Perform a generic Newton-Schulz iteration on a square matrix.

All of the above functions use the iterative approach described below, with several choices of coefficients available.

- `ClassicalNewtonSchulz`: The classical Newton-Schulz coefficients. This method converges slowly compared to other methods, do it is generally not recommended.
- `NSJordan`: The coefficients presented in https://kellerjordan.github.io/posts/muon/. These coefficients are designed to get a fast approximation of `msign` under very low tolerances ($O(0.1)$).
- `PolarExpress` A set of coefficients designed for high tolerances (See https://arxiv.org/abs/2505.16932).
- `NSJianlinSu` A set of coefficients designed in https://kexue.fm/archives/10996 using the same principle as `PolarExpress`.

For general purposes, either `PolarExpress` or `NSJianlinSu` are suitable. The number of steps should in general be tuned to obtain the desired tolerance.


## Matrix Sign and Newton-Schulz

The matrix sign of a rectangular matrix $A \in \mathbb{R}^{m \times n}$ is defined as
$$ \text{msign}(A) = UV^\top = (A A^\top)^{-1/2} A = A (A^\top A)^{-1/2}, $$
where $A = U\Sigma V^\top$ is the (thin) SVD of $A$. In essence $\text{msign}(A)$ is the matrix with the same singular vectors as $A$ but with all non-zero singular values set to one.

Using the SVD decomposition is the simplest method to compute $\text{msign}$, but is computationally expensive, and not efficient on GPUs. An alternative approach consists in computing the following sequence
$$ X_{t+1} = a X_t + b (X_t X_t^\top) X_t + c (X_t X_t^\top)^2 X_t, $$
where $X_0 = A / \|A\|_F$, and $(a,b,c)=(2,-1.5,0.5)$. Each iteration preserves the singular vectors while applying the quintic polynomial $p(x)=ax + bx^3 + cx^5$ to the singular values. It can be shown that $\lim_{n\rightarrow \infty}p^{(n)}(x) = \text{sign}(x)$, and so that $\lim_{t\rightarrow \infty} X_t = \text{msign}(A)$. The advantage of this approach is that each iteration only involves matrix multiplications and additions, so it can be performed efficiently on GPU. 

This package exports functions `msign` and `newton_schulz` to perform these computations.

```julia
A = randn(5, 10)

O_1 = msign(A, MSignSVD()) # compute using SVD

O_2 = msign(A, PolarExpress(), 10) # compute using Newton-Schulz
```

Concretely, `msign` can be computed either via SVD using the `MSignSVD` method, or via Newton-Schulz with any of the coefficients outlined above.

## Alternative for square matrices

For square matrices, we can define an alternative function,
$$ \text{mcsgn}(A) = (A^2)^{-1/2} A = A (A^2)^{-1/2}. $$
This latter function is invariant under similarity transformation ($\text{mcsgn}(PAP^{-1}) = \text{mcsgn}(A)$), and if $A$ is diagonalizable, $\text{mcsgn}(A)$ is the matrix with the same eigenvectors as $A$ and eigenvalues replaced by their sign.

This function can be approximated by a similar iteration as for `msign`, namely
$$ X_{t+1} = a X_t + b X_t^3 + c X_t^5. $$
The same coefficients used for `msign` can be used for this.

```julia
A = rand(10, 10)

O = mcsgn(A, PolarExpress())
```

## Matrix square-root and inverse square-root

The `mcsgn` iteration can be adapted to compute the arithmetic matrix square root of a square positive semidefinite matrix $A$, namely, a square positive semidefinite matrix $B$ such that $B^2 = A$. This can be done in this package via the `msqrt` function

```julia
A = randn(10,10)
A = A * A'
B = msqrt(A, PolarExpress(), 10) # More iterations needed
```

A similar method makes it possible to compute the inverse square root, via the `mrsqrt` function

```julia
C = mrsqrt(A, PolarExpress(), 10)
```

Beware that this method *will* perform poorly when $A$ has eigenvalues close to zero. In many cases, we really just want to compute the product $M A^{-1/2}$, which can be done directly in a slightly stabler way with the `matmul_mrsqrt` function

```julia
A = randn(20,10)
B = randn(10,10)
B = B * B'

C = matmul_mrsqrt(A, B, PolarExpress(), 10)
```

## Applications

TODO

## Acknowledgements

Much of the inspiration for this package came from Jianlin Su's excellent blog posts on the subject, and most of the applications implemented in this package are adapted from there.

- https://kexue.fm/archives/10996
- https://kexue.fm/archives/11025
- https://kexue.fm/archives/11056
- https://kexue.fm/archives/11158

AI generated english translations of these can be found [here](https://datasets.osmarks.net/kexue/site/index.html).

## TODO

- [x] Basic SVD version of `msign`
- [x] Generic Newton-Schulz interface
- [x] `msign` implementation
- [x] `mcsgn` implementation
- [ ] AD Rules for `msign` (following https://kexue.fm/archives/11025, https://kexue-tl.pages.dev/11025-The-Derivative-of-the-msign-Operator)
- [x] Matrix square root and inverse square root
- [ ] Matrix nth root
- [ ] Explicit GPU support