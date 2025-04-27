#=
    Single-body and two-body circuit gates for qubits.
=#

abstract type QubitGate{n} <: AbstractGate end
dim(::QubitGate) = 2
Base.length(::QubitGate{n}) where {n} = n

### Updating the parameters for a gate
export setparams!
"""
    setparams!(gate::QubitGate, params)

Update the parameters that define a QubitGate.
"""
function setparams!(gate::QubitGate, params)
    if !(hasfield(typeof(gate), :params))
        throw(ArgumentError("The gate has no parameters."))
    end
    _setparams!(gate, params)
end
function _setparams!(gate::QubitGate, params)
    gate.params = Float64.(params)
    update!(gate)
end
function update!(gate::QubitGate)
    _update!(gate, gate.gate)
end
_update!(::QubitGate, tensor) = nothing
tensor(gate::QubitGate) = gate.gate

### Constant single qubit gates
export QubitXGate,
    QubitYGate, QubitZGate, QubitIdGate, QubitSGate, QubitSxGate, QubitHadamardGate

"""
    QubitXGate()

Create an X gate that acts on a single qubit.
"""
struct QubitXGate <: QubitGate{1} end
const _qubit_x_gate = SMatrix{2,2,ComplexF64}([0 1; 1 0])
tensor(::QubitXGate) = _qubit_x_gate
QubitNottensor() = QubitXtensor()

"""
    QubitYGate()

Create an Y gate that acts on a single qubit.
"""
struct QubitYGate <: QubitGate{1} end
const _qubit_y_gate = SMatrix{2,2,ComplexF64}([0 -1im; 1im 0])
tensor(::QubitYGate) = _qubit_y_gate

"""
    QubitZGate()

Create an Z gate that acts on a single qubit.
"""
struct QubitZGate <: QubitGate{1} end
const _qubit_z_gate = SMatrix{2,2,ComplexF64}([1 0; 0 -1])
tensor(::QubitZGate) = _qubit_z_gate

"""
    QubitIdGate()

Create an identity gate that acts on a single qubit.
"""
struct QubitIdGate <: QubitGate{1} end
const _qubit_id_gate = SMatrix{2,2,ComplexF64}([1 0; 0 1])
tensor(::QubitIdGate) = _qubit_id_gate

"""
    QubitSGate()

Create an S gate that acts on a single qubit.
"""
struct QubitSGate <: QubitGate{1} end
const _qubit_s_gate = SMatrix{2,2,ComplexF64}([1 0; 0 1im])
tensor(::QubitSGate) = _qubit_s_gate


"""
    QubitSxGate()

Create an Sx gate that acts on a single qubit.
"""
struct QubitSxGate <: QubitGate{1} end
const _qubit_sx_gate = SMatrix{2,2,ComplexF64}(0.5*[1+1im 1-1im; 1-1im 1+1im])
tensor(::QubitSxGate) = _qubit_sx_gate


"""
    QubitHadamardGate()

Create an Hadamard gate that acts on a single qubit.
"""
struct QubitHadamardGate <: QubitGate{1} end
const _qubit_hadamard_gate = SMatrix{2,2,ComplexF64}(sqrt(0.5)*[1 1; 1 -1])
tensor(::QubitHadamardGate) = _qubit_hadamard_gate


### Constant two qubit gates 
export QubitCNOTGate,
    QubitCXGate, QubitCNOTReverseGate, QubitCZGate, QubitSWAPGate, QubitiSWAPGate

"""
    QubitCNOTGate()
    QubitCXGate()

Create a controlled-NOT (controlled-X) gate that acts on two qubits.
The first qubit is the control gate, the second is the target gate.
"""
struct QubitCNOTGate <: QubitGate{2} end
const _qubit_cnot_gate =
    SArray{Tuple{2,2,2,2},ComplexF64}([0 0; 0 1;;; 1 0; 0 0;;;; 1 0; 0 0;;; 0 0; 0 1])
tensor(::QubitCNOTGate) = _qubit_cnot_gate
QubitCXGate() = QubitCNOTGate()

"""
    QubitCNOTReverseGate()

Create a controlled-NOT (controlled-X) gate that acts on two qubits.
The first qubit is the target gate, the second is the control gate.
"""
struct QubitCNOTReverseGate <: QubitGate{2} end
const _qubit_cnot_reverse_gate =
    SArray{Tuple{2,2,2,2},ComplexF64}([0 1; 1 0;;; 0 0; 0 0;;;; 0 0; 0 0;;; 1 0; 0 1])
tensor(::QubitCNOTReverseGate) = _qubit_cnot_reverse_gate

"""
    QubitCZGate()

Create a controlled-Z (controlled-PHASE) gate that acts on two qubits.
"""
struct QubitCZGate <: QubitGate{2} end
const _qubit_cz_gate =
    SArray{Tuple{2,2,2,2},ComplexF64}([1 0; 0 1;;; 0 0; 0 0;;;; 0 0; 0 0;;; -1 0; 0 1])
tensor(::QubitCZGate) = _qubit_cz_gate

"""
    QubitSWAPGate()

Create a SWAP gate that acts on two qubits.
"""
struct QubitSWAPGate <: QubitGate{2} end
const _qubit_swap_gate =
    SArray{Tuple{2,2,2,2},ComplexF64}([1 0; 0 0;;; 0 1; 0 0;;;; 0 0; 1 0;;; 0 0; 0 1])
tensor(::QubitSWAPGate) = _qubit_swap_gate

"""
    QubitiSWAPGate()

Create an iSWAP gate that acts on two qubits.
"""
struct QubitiSWAPGate <: QubitGate{2} end
const _qubit_iswap_gate =
    SArray{Tuple{2,2,2,2},ComplexF64}([1 0; 0 0;;; 0 1im; 0 0;;;; 0 0; 1im 0;;; 0 0; 0 1])
tensor(::QubitiSWAPGate) = _qubit_iswap_gate


### Tunable single qubit gates 
export QubitRxGate, QubitRzGate, QubitRyGate, QubitPhaseGate, QubitRotationGate

# Rx gate
mutable struct QubitRxGate <: QubitGate{1}
    params::Float64
    gate::Array{ComplexF64,2}
end
"""
    QubitRxGate(θ::Number)

Create a rotation gate in the X-axis, exp(-iθX/2),  with angle θ.
"""
function QubitRxGate(θ::Number)
    gate = QubitRxGate(Float64(θ), Array{ComplexF64,2}([0 0; 0 0]))
    update!(gate)
    return gate
end
function _update!(gate::QubitRxGate, tensor)
    tensor[1, 1] = tensor[2, 2] = cos(gate.params/2)
    tensor[1, 2] = tensor[2, 1] = -sin(gate.params/2)*1im
end

# Rz gate
mutable struct QubitRzGate <: QubitGate{1}
    params::Float64
    gate::Array{ComplexF64,2}
end
"""
    QubitRzGate(θ::Number)

Create a rotation gate in the Z-axis, exp(-iθZ/2),  with angle θ.
"""
function QubitRzGate(θ::Number)
    gate = QubitRzGate(Float64(θ), Array{ComplexF64,2}([0 0; 0 0]))
    update!(gate)
    return gate
end
function _update!(gate::QubitRzGate, tensor)
    costerm = cos(gate.params/2)
    sinterm = sin(gate.params/2)*1im
    tensor[1, 1] = costerm - sinterm
    tensor[2, 2] = costerm + sinterm
    tensor[1, 2] = tensor[2, 1] = 0
end

# Ry gate
mutable struct QubitRyGate <: QubitGate{1}
    params::Float64
    gate::Array{ComplexF64,2}
end
"""
    QubitRyGate(θ::Number)

Create a rotation gate in the Y-axis, exp(-iθY/2), with angle θ.
"""
function QubitRyGate(θ::Number)
    gate = QubitRyGate(Float64(θ), Array{ComplexF64,2}([0 0; 0 0]))
    update!(gate)
    return gate
end
function _update!(gate::QubitRyGate, tensor)
    costerm = cos(gate.params/2)
    sinterm = sin(gate.params/2)
    tensor[1, 1] = tensor[2, 2] = costerm
    tensor[1, 2] = -sinterm
    tensor[2, 1] = sinterm
end

# Phase gate
mutable struct QubitPhaseGate <: QubitGate{1}
    params::Float64
    gate::Array{ComplexF64,2}
end
"""
    QubitPhaseGate(θ::Float64)

Create a phase rotation gate with angle θ.
"""
function QubitPhaseGate(θ::Number)
    gate = QubitPhaseGate(Float64(θ), Array{ComplexF64,2}([0 0; 0 0]))
    update!(gate)
    return gate
end
function _update!(gate::QubitPhaseGate, tensor)
    tensor[1, 1] = 1
    tensor[1, 2] = tensor[2, 1] = 0
    tensor[2, 2] = exp(1im*gate.params)
end

mutable struct QubitRotationGate <: QubitGate{1}
    params::NTuple{3,Float64}
    gate::Array{ComplexF64,2}
end
"""
    QubitRotationGate(θ::Number, λ::Number, φ::Number)

Create a general qubit rotation gate (up to arbitary phase). The matrix goes as

```
    [cos(θ/2)           -exp(iλ)sin(θ/2);
     exp(iφ)sin(θ/2)    exp(i(λ+φ))cos(θ/2)]
```
"""
function QubitRotationGate(θ::Number, λ::Number, φ::Number)
    gate = QubitRotationGate(Float64.((θ, λ, φ)), Array{ComplexF64,2}([0 0; 0 0]))
    update!(gate)
    return gate
end
function _update!(gate::QubitRotationGate, tensor)
    costerm = cos(gate.params[1]/2)
    sinterm = sin(gate.params[1]/2)
    tensor[1, 1] = costerm
    tensor[1, 2] = -exp(1im*gate.params[2])*sinterm
    tensor[2, 1] = exp(1im*gate.params[3]) * sinterm
    tensor[2, 2] = exp(1im*(gate.params[2]+gate.params[3])) * costerm
end
