# TeNe.jl

[![Build Status](https://github.com/lcauser/TeNe.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/lcauser/TeNe.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation - DEV](https://img.shields.io/badge/Documentation-DEV-2ea44f)](https://lcauser.github.io/TeNe.jl/dev/)
[![codecov](https://codecov.io/gh/lcauser/TeNe.jl/graph/badge.svg?token=0A2XFLRLMN)](https://codecov.io/gh/lcauser/TeNe.jl)

TeNe.jl is a Julia package for running tensor network calculations.
Please note that it is currently **in development**, but with rapid progress!
A list of some of the features it current includes are
1. Tensor functionality (e.g. tensor contractions, tensor products) with memory caching to reduce the amount of allocations needed.
2. A natural syntax for implementing common and fundamental tensor network operations, e.g., contracting one tensor network with another using syntax such as `ψ * ψ`, or `inner(ψ, ψ)`.
3. High-level interfaces for dealing with physical systems, such as automatic MPO contruction and quantum circuits (in progress).
4. An environment for performing and building variational algorithms for matrix product states. Algorithms such as the density-matrix renormalization group are included by default.

This will remain in development until all the matrix product states functionality is concluded, along with support for CUDA.
At that point we will release v1.0.0, and then add functionality for other tensor network methods, such as uniform matrix product states, tree tensor networks and projected-entangled pair states.


Developing
--------------
If you wish to contribute to TeNe, then your contribution is welcome! TeNe uses JuliaFormatter.jl to ensure it is 
compliant with the Julia guidelines. PRs that are not will be blocked - we use a pre-commit hook to check and 
perform formatting on commits. To use this, please first make sure pre-commit is installed, and run
`pre-commit install`. 