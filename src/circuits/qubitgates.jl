# Qubit Gates 
abstract type QubitGate{n} <: AbstractGate end
dim(::QubitGate) = 2
Base.length(::QubitGate{n}) where {n} = n

export setparams!
function setparams!(gate::QubitGate, params)
    if !(hasfield(typeof(Gate), :params))
        throw(ArgumentError("The gate has no parameters."))
    end
    gate.params = params
end
function _setparams!(gate::QubitGate, params)
    gate.params = params
end


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
export QubitCNOTGate, QubitCXGate, QubitCNOTReverseGate, QubitSWAPGate
struct QubitCNOTGate <: QubitGate{2} end
const _qubit_cnot_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([0 0; 0 1;;; 1 0; 0 0;;;; 1 0; 0 0;;; 0 0; 0 1])
tensor(::QubitCNOTGate) = _qubit_cnot_gate
QubitCXGate() = QubitCNOTGate()

struct QubitCNOTReverseGate <: QubitGate{2} end
const _qubit_cnot_reverse_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([0 1; 1 0;;; 0 0; 0 0;;;; 0 0; 0 0;;; 1 0; 0 1])
tensor(::QubitCNOTReverseGate) = _qubit_cnot_reverse_gate

struct QubitSWAPGate <: QubitGate{2} end
const _qubit_swap_gate = SArray{Tuple{2, 2, 2, 2}, ComplexF64}([1 0; 0 0;;; 0 1; 0 0;;;; 0 0; 1 0;;; 0 0; 0 1])
tensor(::QubitSWAPGate) = _qubit_swap_gate