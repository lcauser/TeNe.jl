#=
    The time-dependant variational principle is used to time-evolve an MPS under
    some Hamiltonian by projecting the evolution into the manifold of MPS.
=#

export tdvp 

"""
    tdvp()
"""
function tdvp()

end

### The updater 
struct MPSUpdateTDVP <: MPSUpdate 
    kryloviters::Int
    krylovdim::Int
    ishermitian::Bool 
    tol::Float64
    dt::Number
end

function MPSUpdateTDVP(
    dt::Number;
    kryloviters::Int=30,
    krylovdim::Int=30,
    krylovtol::Float64=1e-14,
    ishermitian::Bool=true,
    )
    return MPSUpdateTDVP(kryloviters, krylovdim, ishermitian, krylovtol, dt)
end

function update(updater::MPSUpdateTDVP, optim::MPSOptimiser, A)
    # Exponentiate the local site
    f(x) = _tdvp_applyH_M(optim, x)
    Anew = exponentiate(f, -1im*updater.dt, A,
        maxiter=updater.kryloviters,
        krylovdim=updater.krylovdim,
        ishermitian=updater.ishermitian,
        tol=updater.tol
    )

    # Backwards evolution of center 
    if optim.nsites == 1
        
    end

    return Anew
end

function _tdvp_applyH_M(optim::MPSOptimiser, A)
    Anew = similar(A) .= 0
    for projψ in optim.projψs
        Anew .+= product(projψ, A, optim.dir)
    end
    return Anew
end