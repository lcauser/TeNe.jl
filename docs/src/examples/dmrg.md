# Density-matrix renormalization group (DMRG)

The most profound and impactful tensor-network based algorithm is the much celebrated DMRG method,
most-commonly used for estimating the ground-state properties of a one-dimensional quantum Hamiltonian.
The idea is to estimate the ground state by a matrix product state (MPS) $\ket{\psi}$ with bond dimension $\chi$, and update the tensors in the MPS one-by-one until the ground state energy has converged.
In practice, we construct the Hamiltonian $\hat{H}$ as an exact matrix product operator (MPO), and then the energy can be written as 
```math
    E = \frac{\braket{\psi | \hat{H} | \psi}}{\braket{\psi | \psi}}.
```

## Constructing the MPO
Finding the MPO for $\hat{H}$ is not always an easy task.
TeNe is able to construct the MPO for $\hat{H}$ automatically using a simple interface. 
Below, we create the Hamiltonian for a transverse-field Ising model
```math
    \hat{H} = -h\sum_{j=1}^{N} \hat{X}_{i} - J \sum_{j=1}^{N-1} \hat{Z}_{j}\hat{Z}_{j+1}
```
with $h = 1.0$ and $J = 1.0$.
```julia
N = 20
qu = Qubits()
H = OpList(qu, N)
for j = 1:N
    add!(H, "x", j, -1.0)
end
for j = 1:N-1
    add!(H, ["z", "z"], [j, j+1], -1.0)
end
H = MPO(H)
```
The `Qubits()` function creates a `LatticeType` that stores all the key information for qubits, such as different local states (e.g. the "up" and "down" state), and operators (such as "x" or "z").
Then, we create an operator list for $N=20$ qubits, passing through the qubits type so that it knows how to 
identify operators.
We then add the terms to the list like an equation, and convert it to an MPO.

## Finding the ground state 
To find the ground state, we want to call the DMRG algorithm.
We must first have some initial MPS guess, which we'll choose to be random `ψ = randommps(2, N, 1)`, which initalises the MPS with bond dimension $\chi=1$.
This is all we need to call the DMRG algorithm.
```julia
ψ = randommps(2, N, 1)
energy, optim = dmrg(ψ, H; cutoff=1e-12, maxdim=16, maxsweeps=30)
```
There are many keyword arguments you can specify to control the hyperparameters for the method,
and we recommend to look at the manual.
However, we will list some important ones here:
- `cutoff`: The bond dimension of the MPS is controlled by a singular value decomposition with some error tolerance. The cutoff is this tolerance: large values will give a smaller bond dimension (and is thus quicker), but is less precise. On the other hand, a small value gives a high-precision result, but can be slow. Some good values are in the range $10^{-8}$ to $10^{-16}$.
- `maxdim`: The maximum bond dimension for the algorithm. Setting `maxdim=0` means there is no limit on the bond dimension.
- `maxsweeps`: The maximum number of iterations to do.
- `minsweeps`: The minumum number of iterations to run.
A strategy we suggest to use the algorithm with a large `maxdim`, but use the `cutoff` to control the bond dimension, on a scheduele from large-to-small.
```julia
ψ = randommps(2, N, 1)
energy, optim = dmrg(ψ, H; cutoff=1e-6, maxdim=512, maxsweeps=10)
sweep!(optim, 10; cutoff=1e-8)
sweep!(optim, 10; cutoff=1e-10)
sweep!(optim, 10; cutoff=1e-12)
```
Notice that we can continue to do iterations using DMRG using the optimiser `optim`: the `10` is the number of (maximum) number of iterations to do, and the hyperparameters can be changed using key arguments.
The output should look something like :
```juliarepl
iter=1, objective=-2.44896927E+01, maxbonddim=2
iter=2, objective=-2.50832673E+01, maxbonddim=4 
iter=3, objective=-2.51074298E+01, maxbonddim=6
iter=4, objective=-2.51077236E+01, maxbonddim=5
iter=5, objective=-2.51077322E+01, maxbonddim=5
iter=6, objective=-2.51077961E+01, maxbonddim=7
iter=7, objective=-2.51077966E+01, maxbonddim=8
iter=8, objective=-2.51077966E+01, maxbonddim=8
iter=9, objective=-2.51077971E+01, maxbonddim=10
iter=10, objective=-2.51077971E+01, maxbonddim=10
iter=11, objective=-2.51077971E+01, maxbonddim=14 
iter=12, objective=-2.51077971E+01, maxbonddim=14 
```
It is clear that at this point, reducing the `cutoff` any further gives very minor improvement to the ground state energy.
We can see the results below:
```@eval
using TeNe
using Plots

N = 20
qu = Qubits()
H = OpList(qu, N)
for j = 1:N
    add!(H, "x", j, -1.0)
end
for j = 1:N-1
    add!(H, ["z", "z"], [j, j+1], -1.0)
end
H = MPO(H)
ψ = randommps(2, N, 1)
energy, optim = dmrg(ψ, H; cutoff=1e-6, maxdim=512, maxsweeps=10)
sweep!(optim, 10; cutoff=1e-8)
sweep!(optim, 10; cutoff=1e-10)
sweep!(optim, 10; cutoff=1e-12)
plot(optim.costs)
xlabel!("Sweeps")
ylabel!("Energy")
savefig("dmrg_energy.pdf")
nothing
```
![](dmrg_energy.pdf)

## Measuring observables
We now have the ground state and the ground state energy, but these alone are not all that interesting. 
What we would really like to do is measure he properties of the ground state, such as the magnetisations $\hat{X}_{j}$ and $\hat{Z}_{j}\hat{Z}_{j+1}$.
```julia
xs = OpList(qu, N)
for j = 1:N
    add!(xs, "x", j)
end
xs = real(inner(ψ, xs, ψ))

zs = OpList(qu, N)
for j = 1:N-1
    add!(zs, ["z", "z"], [j, j+1])
end
zs = real(inner(ψ, zs, ψ))
```

## Finding excited states
While DMRG is particularly well suited to extremal eigenstates, it can be used reasonably well to estimate the low-lying excited states.
This is done by simply adding the (outer product of the) state `ψ` to our Hamiltonian.
```julia
ψ′ = randommps(2, N, 1)
energy, optim = dmrg(ψ′, H, ψ; cutoff=1e-12, maxdim=128, maxsweeps=30)
```