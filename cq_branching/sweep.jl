# Sweep sigma from 0 to SIGMAX for a fixed mu, one process (no repeated compile).
# Usage: julia -t auto sweep.jl mu sigmax nsteps nsamples seed outdir
include(joinpath(@__DIR__, "sim.jl"))
using Statistics
mu       = parse(Float64, ARGS[1])
sigmax   = parse(Float64, ARGS[2])
nsteps   = parse(Int,     ARGS[3])
nsamples = parse(Int,     ARGS[4])
seed     = parse(Int,     ARGS[5])
outdir   = ARGS[6]
mkpath(outdir)
for (k, sigma) in enumerate(range(0.0, sigmax, length = nsteps))
    t0 = time()
    cq = run(mu, sigma, nsamples, seed + 1000k)
    fn = joinpath(outdir, "cq_mu$(mu)_sig$(lpad(string(round(sigma, digits=4)), 6, '0')).bin")
    open(fn, "w") do io; write(io, cq); end
    println("sigma=$(round(sigma, digits=4))  mean=$(round(mean(cq), digits=4)) std=$(round(std(cq), digits=4))  [$(round(time()-t0, digits=1)) s]")
    flush(stdout)
end
