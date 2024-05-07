using TeNe 

N = 100
J = 1.0
h = 1.0 

qu = Qubits()

X = OpList(qu, N)
for i = 1:N
    add!(X, "x", i, -h)
end
X = MPO(X)

ZZ = OpList(qu, N)
for i = 1:N-1
    add!(ZZ, ["z", "z"], [i, i+1], -J)
end
ZZ = MPO(ZZ)

ψ = randommps(2, N, 1)

energy, optim = dmrg(ψ, X, ZZ; cutoff=1e-10);

println(energy)