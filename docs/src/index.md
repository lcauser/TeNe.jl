# TeNe.jl Documentation

TeNe.jl is a Julia package for running tensor network calculations.
Please note that it is currently **in development**, but with rapid progress!
A list of some of the features it current includes are
1. Tensor functionality (e.g. tensor contractions, tensor products) with memory caching to reduce the amount of allocations needed.
2. A natural syntax for implementing common and fundamental tensor network operations, e.g., contracting one tensor network with another using syntax such as `ψ * ψ`, or `inner(ψ, ψ)`.
3. High-level interfaces for dealing with physical systems, such as automatic MPO contruction and quantum circuits (in progress).
4. An environment for performing and building variational algorithms for matrix product states. Algorithms such as the density-matrix renormalization group are included by default.

This will remain in development until all the matrix product states functionality is concluded, along with support for CUDA.
At that point we will release v1.0.0, and then add functionality for other tensor network methods, such as uniform matrix product states, tree tensor networks and projected-entangled pair states.

This documentation is styled in a way that we hope you can learn from example.
However, we also provide afterwards a "manual" of sorts which details more advanced functionality.
