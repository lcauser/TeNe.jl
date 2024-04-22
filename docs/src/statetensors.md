# State tensors

State tensors are full-rank tensorial objects for describing tensors of a
many-body system, e.g., wavefunctions for lattice systems.

```@contents
Pages = ["statetensors.md"]
Depth = 3
```

## Vectors
### Initiating a state vector

### Properties of state vectors

```@docs
rank(::GStateTensor)
length(::GStateTensor)
norm(::GStateTensor)
```

### Manipulations of state vectors

#### Normalization
```@docs
normalize!(::GStateTensor)
```

### Inner products

### Sampling a state vector

## Operators

### Initiating an operator

### Construct an operator from a list

### Products