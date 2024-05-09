#=
    Single-body and two-body circuit gates for qubits.
=#

abstract type QubitGate{n} <: AbstractGate end
dim(::QubitGate) = 2
Base.length(::QubitGate{n}) where {n} = n

### Updating the parameters for a gate
export setparams!
function setparams!(gate::QubitGate, params)
    if !(hasfield(typeof(Gate), :params))
        throw(ArgumentError("The gate has no parameters."))
    end
    _setparams!(gate, params)
end
function _setparams!(gate::QubitGate, params)
    gate.params = params
    update!(gate)
end
function update!(gate::QubitGate)
    _update!(gate, gate.gate)
end
_update!(gate, tensor) = nothing
tensor(gate::QubitGate) = gate.gate

### Constant single qubit gates
export QubitXGate, QubitYGate, QubitZGate, QubitIdGate, QubitSGate, QubitSxGate, QubitHadamardGate

struct QubitXGate <: QubitGate{1} end
const _qubit_x_gate = SMatrix{2, 2, ComplexF64}([0 1; 1 0])
tensor(::QubitXGate) = _qubit_x_gate
QubitNottensor() = QubitXtensor()

struct QubitYGate <: QubitGate{1} end
const _qubit_y_gate = SMatrix{2, 2, ComplexF64}([0 -1im; 1im 0])
tensor(::QubitYGate) = _qubit_y_gate

struct QubitZGate <: QubitGate{1} end
const _qubit_z_gate = SMatrix{2, 2, ComplexF64}([1 0; 0 -1])
tensor(::QubitZGate) = _qubit_z_gate

struct QubitIdGate <: QubitGate{1} end
const _qubit_id_gate = SMatrix{2, 2, ComplexF64}([1 0; 0 1])
tensor(::QubitIdGate) = _qubit_id_gate

struct QubitSGate <: QubitGate{1} end
const _qubit_s_gate = SMatrix{2, 2, ComplexF64}([1 0; 0 1im])
tensor(::QubitSGate) = _qubit_s_gate

struct QubitSXGate <: QubitGate{1} end
const _qubit_sx_gate = SMatrix{2, 2, ComplexF64}(0.5*[1+1im 1-1im; 1-1im 1+1im])
tensor(::QubitSXGate) = _qubit_sx_gate

struct QubitHadamardGate <: QubitGate{1} end
const _qubit_hadamard_gate = SMatrix{2, 2, ComplexF64}(sqrt(0.5)*[1 1; 1 -1])
tensor(::QubitHadamardGate) = _qubit_hadamard_gate


### Constant two qubit gates 
export QubitCNOTGate, QubitCXGate, QubitCNOTReverseGate, QubitSWAPGate, QubitiSWAPGate
struct QubitCNOTGate <: QubitGate{2} end
const _qubit_cnot_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([0 0; 0 1;;; 1 0; 0 0;;;; 1 0; 0 0;;; 0 0; 0 1])
tensor(::QubitCNOTGate) = _qubit_cnot_gate
QubitCXGate() = QubitCNOTGate()

struct QubitCNOTReverseGate <: QubitGate{2} end
const _qubit_cnot_reverse_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([0 1; 1 0;;; 0 0; 0 0;;;; 0 0; 0 0;;; 1 0; 0 1])
tensor(::QubitCNOTReverseGate) = _qubit_cnot_reverse_gate

struct QubitCZGate <: QubitGate{2} end
const _qubit_cz_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([1 0; 0 1;;; 0 0; 0 0;;;; 0 0; 0 0;;; -1 0; 0 1])
tensor(::QubitCZGate) = _qubit_cz_gate 

struct QubitSWAPGate <: QubitGate{2} end
const _qubit_swap_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([1 0; 0 0;;; 0 1; 0 0;;;; 0 0; 1 0;;; 0 0; 0 1])
tensor(::QubitSWAPGate) = _qubit_swap_gate

struct QubitiSWAPGate <: QubitGate{2} end 
const _qubit_iswap_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([1 0; 0 0;;; 0 1im; 0 0;;;; 0 0; 1im 0;;; 0 0; 0 1])
tensor(::QubitiSWAPGate) = _qubit_swap_gate


### Tunable single qubit gates 
export QubitRxGate, QubitRzGate, QubitRyGate, QubitPhaseGate

# Rx gate
struct QubitRxGate{Q} <: QubitGate{1} where {Q<:Number}
    params::Q
    gate::SArray{Tuple{2, 2}, ComplexF64}
end
function _update!(gate::QubitRxGate, tensor)
    tensor[1, 1] = tensor[2, 2] = cos(gate.params/2)
    tensor[1, 2] = tensor[2, 1] = -sin(gate.params/2)*1im
end

# Rz gate
struct QubitRzGate{Q} <: QubitGate{1} where {Q<:Number}
    params::Q
    gate::SArray{Tuple{2, 2}, ComplexF64}
end
function _update!(gate::QubitRzGate, tensor)
    costerm = cos(gate.params/2)
    sinterm = sin(gate.params/2)*1im
    tensor[1, 1] = costerm - sinterm
    tensor[2, 2] = costerm + sinterm
    tensor[1, 2] = tensor[2, 1] = 0
end

# Ry gate
struct QubitRyGate{Q} <: QubitGate{1} where {Q<:Number}
    params::Q
    gate::SArray{Tuple{2, 2}, ComplexF64}
end
function _update!(gate::QubitRyGate, tensor)
    costerm = cos(gate.params/2)
    sinterm = sin(gate.params/2)
    tensor[1, 1] = tensor[2, 2] = costerm 
    tensor[1, 2] = -sinterm
    tensor[2, 1] = sinterm
end

# Phase gate
struct QubitPhaseGate{Q} <: QubitGate{1} where {Q<:Number}
    params::NTuple{3, Q}
    gate::SArray{Tuple{2, 2}, ComplexF64}
end
function _update!(gate::QubitPhaseGate, tensor)
    tensor[1, 1] = 1
    tensor[1, 2] = tensor[2, 1] = 0
    tensor[2, 2] = exp(1im*gate.params)
end

struct QubitRotationGate{Q} where {Q}
    params::NTuple{3, Q}
    gate::SArray{Tuple{2, 2}, ComplexF64}
end