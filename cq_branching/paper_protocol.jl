# Exact re-implementation of the authors' MATLAB protocol
# (github.com/Azuresky99/quPCR: Fig1D_SimSingleDnaPCR.m, FigS3B_SimSingleDnaRandomE.m).
#
#   * Stochastic stage: 19 binomial doublings (their loop roundNum = 2:20).
#     Per tube and per cycle one efficiency e_n ~ N(1+AE, s) is drawn
#     (shared by all molecules in that tube); a molecule replicates if
#     rand <= e_n - 1, so p is effectively clipped to [0,1].
#   * Deterministic stage: Cq = log(T / N_20) / log(1+AE) + 19,
#     with T = (1+AE)^quCq. No noise after cycle 20.
#   * FigS3B: AE = 0.9593, quCq = 38.335, s = 0.1 / z_{0.975} = 0.0510
#     ("± 0.1" is the half-width of a 95% interval, not one sigma).
#
# Usage: julia -t auto paper_protocol.jl AE quCq halfwidth nstoch nsamples seed out.bin [noise_all_cycles]
#   noise_all_cycles = 1 also applies the per-cycle noise after cycle 20
#   (log-growth random walk, same s), which the authors' code omits.
include(joinpath(@__DIR__, "sim.jl"))
using Random, Statistics

function paper_cq(rng, ae, qucq, s, nstoch, noise_late)
    N = 1
    for _ in 1:nstoch
        p = s == 0.0 ? ae : (1.0 + ae + s * randn(rng)) - 1.0
        N += binom(rng, N, clamp(p, 0.0, 1.0))
    end
    logT = qucq * log1p(ae)
    cq = (logT - log(N)) / log1p(ae) + nstoch
    if noise_late && s > 0
        # remaining cycles: log N grows by log(1+p_n); find first passage of log T
        x = log(N); n = nstoch
        while true
            n += 1
            g = log1p(clamp(ae + s * randn(rng), 0.0, 1.0))
            if x + g >= logT
                return (n - 1) + (logT - x) / g
            end
            x += g
        end
    end
    return cq
end

if abspath(PROGRAM_FILE) == @__FILE__
    ae, qucq, hw = parse.(Float64, ARGS[1:3])
    nstoch, nsamp, seed = parse.(Int, ARGS[4:6])
    out = ARGS[7]
    late = length(ARGS) >= 8 && ARGS[8] == "1"
    s = hw / 1.959963984540054
    cq = Vector{Float64}(undef, nsamp)
    rngs = [Xoshiro(seed + 7919t) for t in 1:Threads.nthreads()]
    Threads.@threads :static for i in 1:nsamp
        cq[i] = paper_cq(rngs[Threads.threadid()], ae, qucq, s, nstoch, late)
    end
    open(out, "w") do io; write(io, cq); end
    println("AE=$ae quCq=$qucq sigma=$(round(s, digits=4)) late_noise=$late  mean=$(round(mean(cq), digits=3)) std=$(round(std(cq), digits=3))")
end
