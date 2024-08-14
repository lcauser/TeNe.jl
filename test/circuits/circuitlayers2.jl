N = 20
J = 1.0
h = 1.0
# Create Hamiltonian 
qu = Qubits()
H = OpList(qu, N)
for i = 1:N
    add!(H, "x", i, -h)
end
for i = 1:N-1
    add!(H, ["z", "z"], [i, i+1], -J)
end
H = MPO(H)

# Do DMRG 
ψ = randommps(2, N, 1)
energy, optim = dmrg(ψ, H; cutoff=1e-12, verbose=true)


# Create circuit 
depth = 4
circuit = randombwcircuit(2, N, depth)
ϕ = productmps(qu, ["up" for _ = 1:N])
projU = ProjMPSCircuit(ψ, circuit, ϕ);

m = 1
g = 1

# Move the center 
movecenter!(projU, depth+1-m)

# Get the MPS evolved 
top = deepcopy(TeNe.topblock(projU, depth-m))
applygates!(top, projU.circuit.layers[m]; cutoff=1e-16)
gate = deepcopy(projU.circuit.layers[m].gates[g])

# Undo one of the gates
qubits = projU.circuit.layers[m].sites[g]
applygate!(gate, top, qubits[1]; cutoff=1e-16)

# Get the bottom and make a projection
bottom = TeNe.bottomblock(projU, depth-m+2)
projmps = ProjMPS(top, bottom; center=qubits[1])