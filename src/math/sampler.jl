struct UniformSampler
    rng::ChaChaStream

    UniformSampler() = new(ChaCha12Stream())
end

"""
uniform_binary(N, hw = h) outputs a uniform binary Int64 vector with Hamming weight h.
If hw is not given or set to zero, it outputs a truly uniform binary vector.
"""
uniform_binary(us::UniformSampler, N::Int64, hw::Int64=0) = begin
    if hw == 0
        rand(us.rng, [0, 1], N)
    else
        res = zeros(Int64, N)
        cnt = 0
        while cnt < hw
            idx = rand(us.rng, 1:N)
            if res[idx] == 0
                res[idx] = 1
                cnt += 1
            end
        end
        res
    end
end

"""
uniform_ternary(N, hw = h) outputs a uniform ternary Int64 vector with Hamming weight h.
If hw is not given or set to zero, it outputs a truly uniform ternary vector.
"""
uniform_ternary(us::UniformSampler, N::Int64, hw::Int64=0) = begin
    @assert 0 ≤ hw ≤ N "The Hamming weight must be between 0 and N."

    if hw == 0
        rand(us.rng, [-1, 0, 1], N)
    else
        res = zeros(Int64, N)
        cnt = 0
        while cnt < hw
            idx = rand(us.rng, 1:N)
            if res[idx] == 0
                res[idx] = rand(us.rng, -1:2:1)
                cnt += 1
            end
        end
        res
    end
end

"""
uniform_random outputs a uniform random distribution.
"""
uniform_random(us::UniformSampler, N::Int64, Q::Modulus) = rand(us.rng, 0:Q.Q-1, N)

uniform_random(us::UniformSampler, N::Int64, Q::Moduli) = begin
    res = Vector{Vector{UInt64}}(undef, length(Q))
    @inbounds for i = eachindex(Q)
        res[i] = Vector{UInt64}(undef, N)
        for j = 1:N
            res[i][j] = rand(us.rng, 0:Q[i].Q-1)
        end
    end

    res
end

uniform_random(us::UniformSampler, Q::Modulus) = rand(us.rng, 0:Q.Q-1)

uniform_random_to!(us::UniformSampler, x::AbstractVector{UInt64}, Q::Modulus) = begin
    N = length(x)
    @inbounds for i = 0:8:N-8
        x[i+1] = rand(us.rng, 0:Q.Q-1)
        x[i+2] = rand(us.rng, 0:Q.Q-1)
        x[i+3] = rand(us.rng, 0:Q.Q-1)
        x[i+4] = rand(us.rng, 0:Q.Q-1)
        x[i+5] = rand(us.rng, 0:Q.Q-1)
        x[i+6] = rand(us.rng, 0:Q.Q-1)
        x[i+7] = rand(us.rng, 0:Q.Q-1)
        x[i+8] = rand(us.rng, 0:Q.Q-1)
    end
    @inbounds for i = 8(N>>3)+1:N
        x[i] = rand(us.rng, 0:Q.Q-1)
    end
end

uniform_random_to!(us::UniformSampler, x::AbstractVector{Vector{UInt64}}, Q::Moduli) = begin
    @assert length(x) == length(Q) "The array size does not match the moduli."
    @inbounds for i = eachindex(Q)
        N = length(x[i])
        for j = 0:8:N-8
            x[i][j+1] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+2] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+3] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+4] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+5] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+6] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+7] = rand(us.rng, 0:Q[i].Q-1)
            x[i][j+8] = rand(us.rng, 0:Q[i].Q-1)
        end
        for j = 8(N>>3)+1:N
            x[i][j] = rand(us.rng, 0:Q[i].Q-1)
        end
    end
end

export UniformSampler, uniform_binary, uniform_ternary, uniform_random, uniform_random_to!

#=============================================================================================#

samplebit(rng::ChaChaStream) = rand(rng, UInt8) & one(UInt64)

const tailcut = 6

struct CDTSampler
    rng::ChaChaStream       # Base sampler.
    table::Vector{UInt64}   # Probability Table.
    center::Float64         # Distribution center
    sigma::Float64          # Standard Deviation
    taillow::Int64
    tailhigh::Int64
    cint::Int64
    cfrac::Float64

    function CDTSampler(center::Real, sigma::Real)
        normFactor = typemax(UInt64)
        cInt = floor(Int64, center)
        cFrac = center - cInt

        tailLow = round(Int64, center - tailcut * sigma)
        tailHigh = round(Int64, center + tailcut * sigma)
        tailCount = tailHigh - tailLow + 1

        table = Vector{UInt64}(undef, tailCount)
        cdf = 0.0

        for i = 0:tailHigh-tailLow
            xf = Float64(i + tailLow)
            rho = exp(-π * (xf - cFrac)^2 / sigma^2) / sigma
            cdf += rho
            table[i+1] = cdf ≥ 1 ? typemax(UInt64) : round(UInt64, cdf * normFactor)
        end

        new(ChaCha12Stream(), table, center, sigma, tailLow, tailHigh, cInt, cFrac)
    end
end

sample(s::CDTSampler) = searchsortedlast(s.table, rand(s.rng, UInt64)) + s.cint + s.taillow

const baselog = 10
const sampledepth = 3
const hipreclog = baselog * sampledepth
const lowpreclog = 53 - hipreclog

struct VarCenSampler
    basesamplers::Vector{CDTSampler}

    function VarCenSampler(sigma::Real)
        bk = 0.0
        for i = 0:sampledepth-1
            bk += Float64(1 << baselog)^(-2i)
        end
        sigma /= bk

        basesamplers = Vector{CDTSampler}(undef, 1 << baselog)
        for i = 0:1<<baselog-1
            c = i / (1 << baselog)
            basesamplers[i+1] = CDTSampler(c, sigma)
        end

        new(basesamplers)
    end
end

function sample(s::VarCenSampler, center::Float64)
    cInt = floor(Int64, center)
    cFrac = center - cInt
    cFrac64 = unsafe_trunc(UInt64, cFrac * (1 << 53))

    cFrac64Hi = (cFrac64 >> lowpreclog) % Int64
    r = rand(s.basesamplers[1].rng, UInt64)
    @inbounds for i = reverse(0:lowpreclog-1)
        b = (r >> i) & 1
        cFracBit = (cFrac64 >> i) & 1
        b > cFracBit && (return sampleC(s, cFrac64Hi) + cInt)
        b < cFracBit && (return sampleC(s, cFrac64Hi + 1) + cInt)
    end

    sampleC(s, cFrac64Hi + 1) + cInt
end

function sampleC(s::VarCenSampler, c::Int64)
    mask = (1 << baselog) - 1

    @inbounds for i = 0:sampledepth-1
        r = sample(s.basesamplers[c&mask+1])
        (c & mask > 0 && c < 0) && (r -= 1)
        c >>= baselog
        c += r
    end

    c
end

const cdtbase = 256

struct TwinCDTSampler
    rng::ChaChaStream
    table::Vector{Vector{UInt64}}
    tail_lo::Int64
    tail_hi::Int64
    σ::Float64

    function TwinCDTSampler(σ::Real)
        tables = Vector{Vector{UInt64}}(undef, cdtbase)
        for i = 0:cdtbase-1
            tables[i+1] = compute_cdt(i / cdtbase, σ)
        end

        tail_hi = ceil(Int64, tailcut * σ)
        tail_lo = -tail_hi

        new(ChaCha12Stream(), tables, tail_lo, tail_hi, σ)
    end
end

function compute_cdt(c::Real, σ::Real)
    tail_hi = ceil(Int64, tailcut * σ)
    tail_lo = -tail_hi
    tablesize = tail_hi - tail_lo + 1

    table = Vector{UInt64}(undef, tablesize)
    cdf = 0.0
    for i = 0:tablesize-1
        x = Float64(i + tail_lo)
        rho = exp(-π * (x - c)^2 / σ^2) / σ
        cdf += rho
        if cdf ≥ 1
            table[i+1] = typemax(UInt64)
        else
            table[i+1] = round(UInt64, cdf * typemax(UInt64))
        end
    end

    table
end

function sample(centre::Float64, s::TwinCDTSampler)
    cfloor = floor(centre)
    c0 = floor(Int64, cdtbase * (centre - cfloor))
    c1 = ceil(Int64, cdtbase * (centre - cfloor))

    u = rand(s.rng, UInt64)
    v0 = searchsortedlast(s.table[c0%cdtbase+1], u)
    v1 = searchsortedlast(s.table[c1%cdtbase+1], u)

    if v0 == v1
        return unsafe_trunc(Int64, v0 + cfloor + s.tail_lo)
    end

    cfrac = centre - cfloor
    cdf = 0.0
    for x = s.tail_lo:v0
        xf = Float64(x)
        cdf += exp(-π * (xf - cfrac)^2 / s.σ^2) / s.σ
    end

    p = u / typemax(UInt64)
    if p < cdf
        return unsafe_trunc(Int64, v0 + cfloor + s.tail_lo)
    end
    return unsafe_trunc(Int64, v1 + cfloor + s.tail_lo)
end

struct RGSampler
    rng::ChaChaStream
    τ::Float64

    function RGSampler()
        new(ChaCha12Stream(), 3.2)
    end

    function RGSampler(τ::Real)
        new(ChaCha12Stream(), τ)
    end
end

sample(s::RGSampler) = round(Int64, randn(s.rng) * s.τ / √(2π))

sample(s::RGSampler, τ::Real) = round(Int64, randn(s.rng) * τ / √(2π))

sample128(s::RGSampler, τ::Real) = round(Int128, randn(s.rng) * τ / √(2π))

export CDTSampler, VarCenSampler, TwinCDTSampler, RGSampler, sample