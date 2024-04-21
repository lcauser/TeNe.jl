# Matrix product states

Matrix product states are an ansatz for describing and estimating states of one-dimensional systems.
In this package, MPS can be used for a variety of applications, e.g., ground state estimation, simulating dynamics, and quantum circuit calculations. 

```@contents
Pages = ["mps.md"]
Depth = 3
```

## Matrix product states (MPS)
### Initiating an MPS
There are many ways to create an MPS, and the method will depend on the problem at hand.

#### Random states
For variational problems, we sometimes recommend initiating an MPS randomly using `ψ = randommps(dim, length, bonddim)`, where `dim` is the physical dimension of the lattice, `length` is the number of sites in the lattice, and `bonddim` is the bond dimension of the MPS.

```@docs
randommps
```

#### Product states
An alternative option is to initalise the MPS is a state with known desirable properties. 
This can be done in a few ways.
The simplest way is to just provide a translationally invariant tensor `A` for the MPS (with bond dimension one), `ψ = productmps(N, A)`.

```@docs
productmps(N::Int, A::Q) where {Q<:AbstractArray}
```

### Properties of an MPS
```@docs
rank(::GMPS)
dim(::GMPS)
center(::GMPS)
bonddim(::GMPS, ::Int)
maxbonddim(::GMPS)
norm(::GMPS)
```

### Manipulations of an MPS

#### Normalization
```@docs
normalize!(::GMPS)
```

#### Canonical center
An important property of MPS is that they can be brought into a canonical representation, allowing for better conditioned calculations and simplifications in many MPS algorithms.
The canonical center of an MPS `ψ` is easily moved to `idx` using `movecenter!(ψ, idx)`.
This also allows for dynamic truncation of the MPS, using keyword arguments such as `maxdim` or `cutoff`.
```@docs
movecenter!(ψ::GMPS, idx::Int)
```

#### Truncations
The bond dimension of the MPS controls its accuracy. A larger bond dimension has a higher capacity for precision, but also increases the computational cost of algorithms.
To reduce the bond dimension (and sacrificing as little accuracy as possible), one can truncate the MPS.

```@docs
truncate!(::GMPS)
```

### Inner products
Some linear algebra operations such as the inner product are easy to calculate for MPS.
```@docs
inner(ψ::MPS, ϕ::MPS)
```

## Matrix product operators (MPO)
Just as a wavefunction, or state vector, can be represented as an MPS, an operator can be represented by a matrix product operator.

### Initiating an MPO 


### Construct an MPO from an operator list

### Manipulate an MPO

### Inner products

## Advanced usage: Generalised matrix product states (GMPS)