```@meta 
CollapsedDocStrings = true
```
# State tensors

State tensors are full-rank tensorial objects for describing tensors of a
many-body system, e.g., wavefunctions for lattice systems.

```@contents
Pages = ["statetensors.md"]
Depth = 3
```

## Vectors
### Initiating a state vector
There are many ways to initialise a state. For variational methods, a popular choice is to initiate it randomly, which can be achieved using `randomsv(dim, length)` for physical dimension `dim` and `length` lattice sites.
```@docs
randomsv(::Int, ::Int)
```

Alternatively, you can initalise it by a product state.
```@docs
productsv
```

### Properties of state vectors

```@docs
rank(::GStateTensor)
dim(::GStateTensor)
length(::GStateTensor)
norm(::GStateTensor)
entropy(::StateVector, ::Int)
entropy(::StateVector, ::Any)
```

### Manipulations of state vectors

#### Normalization
```@docs
normalize!(::GStateTensor)
```

#### Exponentiation 
```@docs
exp(O::StateOperator)
```

### Inner products
```@docs
inner(::StateVector, ::StateVector)
```

### Sampling a state vector
```@docs
sample(::StateVector)
```

## Operators
Operators, such as the Hamiltonian of a quantum many-body system, or more generally a matrix can be represented by a StateOperator.

### Initiating an operator
Like a StateVector, we can initalise a StateOperator randomly or as a product state.
```@docs
randomso
productso
```

### Construct an operator from a list

### Products
A StateVector or a StateOperator can be multiplied by a StateOperator.
```@docs
applyso
```

This can be done in-place, either so store the result in an existing StateVector `ϕ`, or to replace `ψ`. In the case of a StateOperator-StateOperator multiplication, the StateOperator for the result must be specified.
```@docs
applyso!
```

The expectation value of operators (or a string of operators) can be calculated in the following way.
```@docs 
inner(::StateVector, ::StateOperator, ::StateVector)
```

The trace of an operator (or a string of operators) can be computed in the following way.
```@docs 
trace(::StateOperator...)
```