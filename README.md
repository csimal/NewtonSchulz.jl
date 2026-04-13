# NewtonSchulz

[![Build Status](https://github.com/csimal/NewtonSchulz.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/csimal/NewtonSchulz.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/csimal/NewtonSchulz.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/csimal/NewtonSchulz.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

A Julia implementation of the Newton-Schulz and its variants for computing the matrix sign function, as well as matrix square roots.

## TODO

- [x] Basic SVD version of `msign`
- [x] Generic Newton-Schulz interface
- [x] `msign` implementation
- [x] `mcsgn` implementation
- [ ] AD Rules for `msign` (following https://kexue.fm/archives/11025, https://kexue-tl.pages.dev/11025-The-Derivative-of-the-msign-Operator)
- [x] Matrix square root and inverse square root
- [ ] Matrix nth root
- [ ] Explicit GPU support