# Stochastic PCR / branching-process simulation of the quantification cycle Cq.
#
# Model: N_0 = 1. At cycle n each particle independently duplicates with
# probability p_n = mu + sigma * xi_n, xi_n ~ N(0,1), drawn once per cycle
# (sigma = 0 reproduces the fixed-AE case).  N_n = N_{n-1} + Binomial(N_{n-1}, p_n).
# A Gaussian p_n can leave [0,1]; two conventions are implemented (PMODE):
#   "clip"   : p_n -> clamp(p_n, 0, 1)                          (default)
#   "extend" : p_n > 1 means each particle makes floor(p_n) copies plus one
#              more with probability p_n - floor(p_n); p_n < 0 -> 0.
#              E[N_n | N_{n-1}] = (1 + p_n) N_{n-1} is then unbiased.
# Stop at the first cycle n with N_n >= T = 2^38.
# Continuous Cq: log-linear interpolation between the last two cycles,
#   Cq = (n-1) + log(T / N_{n-1}) / log(N_n / N_{n-1}).
#
# Only the standard library is used (the package server is unreachable here),
# so the binomial sampler is hand-written:
#   * N <= NEXACT: exact, by geometric jumps over the rarer outcome
#     (cost O(N * min(p, 1-p)) instead of O(N));
#   * N  > NEXACT: Gaussian approximation (relative fluctuation ~ N^{-1/2}
#     is far below anything visible in Cq at that size).
#
# Usage:  julia -t auto sim.jl mu sigma nsamples seed outfile.bin [clip|extend]
# Output: little-endian Float64 array of Cq values.

using Random

const T      = 2.0^38
const NEXACT = 1 << 16
const PMODE  = Ref(:clip)

"Exact Binomial(N, p) via geometric jumps over the rarer outcome."
@inline function binom_exact(rng::AbstractRNG, N::Int, p::Float64)
    p <= 0.0 && return 0
    p >= 1.0 && return N
    # count the rarer outcome, q = its probability
    count_success = p <= 0.5
    q = count_success ? p : 1.0 - p
    lq = log1p(-q)                      # log(1-q) < 0
    k = 0
    i = 0
    while true
        # gap to the next rare event: Geometric(q) on {1,2,...}
        i += ceil(Int, log(rand(rng)) / lq)
        i > N && break
        k += 1
    end
    return count_success ? k : N - k
end

"Binomial(N, p): exact for N <= NEXACT, Gaussian approximation above."
@inline function binom(rng::AbstractRNG, N::Int, p::Float64)
    if N <= NEXACT
        return binom_exact(rng, N, p)
    else
        p <= 0.0 && return 0
        p >= 1.0 && return N
        m = N * p
        s = sqrt(N * p * (1.0 - p))
        k = round(Int, m + s * randn(rng))
        return clamp(k, 0, N)
    end
end

"One trajectory; returns the interpolated Cq."
function one_cq(rng::AbstractRNG, mu::Float64, sigma::Float64)
    N = 1
    n = 0
    while true
        n += 1
        p = sigma == 0.0 ? mu : mu + sigma * randn(rng)
        if PMODE[] === :clip || p <= 1.0
            Nnew = N + binom(rng, N, clamp(p, 0.0, 1.0))
        else                                  # :extend, p > 1
            m = floor(p)
            Nnew = N + Int(m) * N + binom(rng, N, p - m)
        end
        if Nnew >= T
            Nn   = Float64(N)
            Nnn  = Float64(Nnew)
            return (n - 1) + log(T / Nn) / log(Nnn / Nn)
        end
        N = Nnew
        n > 100_000 && return NaN   # p stuck at 0: give up (never happens for the cases run here)
    end
end

function run(mu, sigma, nsamples, seed)
    cq = Vector{Float64}(undef, nsamples)
    nth = Threads.nthreads()
    rngs = [Xoshiro(seed + 7919 * t) for t in 1:nth]
    Threads.@threads :static for i in 1:nsamples
        rng = rngs[Threads.threadid()]
        cq[i] = one_cq(rng, mu, sigma)
    end
    return cq
end

if abspath(PROGRAM_FILE) == @__FILE__
    mu       = parse(Float64, ARGS[1])
    sigma    = parse(Float64, ARGS[2])
    nsamples = parse(Int,     ARGS[3])
    seed     = parse(Int,     ARGS[4])
    outfile  = ARGS[5]
    length(ARGS) >= 6 && (PMODE[] = Symbol(ARGS[6]))
    t0 = time()
    cq = run(mu, sigma, nsamples, seed)
    open(outfile, "w") do io
        write(io, cq)
    end
    using Statistics
    println("mu=$mu sigma=$sigma mode=$(PMODE[]) n=$nsamples threads=$(Threads.nthreads())  ",
            "mean=$(round(mean(cq), digits=4)) std=$(round(std(cq), digits=4)) ",
            "min=$(round(minimum(cq), digits=3)) max=$(round(maximum(cq), digits=3))  ",
            "[$(round(time()-t0, digits=1)) s]")
end
