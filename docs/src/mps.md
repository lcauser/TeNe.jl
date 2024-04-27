```@meta 
CollapsedDocStrings = true
```
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

Alternatively, we can use the LatticeTypes interface to write a product state.
```@docs
productmps(lt::LatticeTypes, states::AbstractVector{String})
```

#### From StateVectors
If you need to write an MPS using some known StateVector, you can call `MPS(ψ::StateVector)`.
```@docs
MPS(::GStateTensor{1})
```

### Properties of an MPS
```@docs
TeNe.rank(::GMPS)
dim(::GMPS)
center(::GMPS)
bonddim(::GMPS, ::Int)
maxbonddim(::GMPS)
norm(::GMPS)
entropy(::MPS, ::Int)
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

### Sampling an MPS
By considering an MPS to be a wavefunction from quantum mechanics (we sometimes call this a *Born machine* in a classical context), we are able to sample it exactly.
```@docs
sample(::MPS)
``` 

## Matrix product operators (MPO)
Just as a wavefunction, or state vector, can be represented as an MPS, an operator can be represented by a matrix product operator.

### Initiating an MPO 
Like the MPS, an MPO can be initiated randomly or as a product operator.
It can also be initiated from some StateOperator.
```@docs
randommpo
productmpo(N::Int, A::AbstractArray; T::Type=ComplexF64)
MPO(::GStateTensor{2})
```

Alternatively, we can use the LatticeTypes interface to write an operator.
```@docs
productmpo(lt::LatticeTypes, ops::AbstractVector{String})
```

### Construct an MPO from an operator list

### Products
The expectation value of a string of MPOs with respect to some MPSs can be calculated exactly. The below can be done with many MPOs.
```@docs
inner(ψ::MPS, O::MPO, ϕ::MPS)
```

Similarly, the trace of a string of MPOs can be calculated.
```@docs
trace(Os::MPO...)
```
## Advanced usage: Generalised matrix product states (GMPS)