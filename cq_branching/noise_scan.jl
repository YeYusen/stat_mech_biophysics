# Noise scan in the authors' protocol for a given AE, at the SAME physical threshold
# T = 1.959^38.34 copies (so quCq(AE) = 38.34 ln 1.959 / ln(1+AE)).
# Usage: julia -t auto noise_scan.jl AE sigmax dsig nsamples outdir
include(joinpath(@__DIR__, "paper_protocol.jl"))
ae = parse(Float64, ARGS[1]); sigmax = parse(Float64, ARGS[2]); dsig = parse(Float64, ARGS[3])
nsamp = parse(Int, ARGS[4]); outdir = ARGS[5]; mkpath(outdir)
qucq = 38.34 * log(1.959) / log1p(ae)
println("AE=$ae quCq=$(round(qucq, digits=3))")
for late in (false, true), sig in 0.0:dsig:sigmax
    fn = joinpath(outdir, "scan_ae$(ae)_late$(Int(late))_sig$(round(sig, digits=4)).bin")
    isfile(fn) && continue                      # resume: skip finished points
    cq = Vector{Float64}(undef, nsamp)
    rngs = [Xoshiro(29 + 7919t + round(Int, 1e5sig) + 13late + round(Int, 1e3ae)) for t in 1:Threads.nthreads()]
    Threads.@threads :static for i in 1:nsamp
        cq[i] = paper_cq(rngs[Threads.threadid()], ae, qucq, sig, 19, late)
    end
    open(fn * ".tmp", "w") do io; write(io, cq); end
    mv(fn * ".tmp", fn; force = true)
end
