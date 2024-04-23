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

### Inner products
```@docs
inner(::StateVector, ::StateVector)
```

### Sampling a state vector
```@docs
sample(::StateVector)
```

## Operators

### Initiating an operator

### Construct an operator from a list

### Products