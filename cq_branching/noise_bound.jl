# Scan the per-cycle AE noise in the authors' protocol and record the Cq histograms,
# to find the largest noise for which peaks 3 and 4 survive.
# Usage: julia -t auto noise_bound.jl nsamples outdir
include(joinpath(@__DIR__, "paper_protocol.jl"))
nsamp = parse(Int, ARGS[1]); outdir = ARGS[2]; mkpath(outdir)
for late in (false, true), sig in 0.0:0.0025:0.05
    cq = Vector{Float64}(undef, nsamp)
    rngs = [Xoshiro(17 + 7919t + round(Int, 1e5sig) + 13late) for t in 1:Threads.nthreads()]
    Threads.@threads :static for i in 1:nsamp
        cq[i] = paper_cq(rngs[Threads.threadid()], 0.9593, 38.335, sig, 19, late)
    end
    open(joinpath(outdir, "scan_late$(Int(late))_sig$(round(sig, digits=4)).bin"), "w") do io; write(io, cq); end
end
