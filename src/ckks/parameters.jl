"""
    CKKSParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, scaling_factor::Integer; dlen::Int64=0)

CKKSParamSketch is a struct for the parameters for the CKKS scheme.
The following parameters stand for:
- `ring_param`: base ring parameters. Currently, only power-of-two cyclotomic rings are supported.
- `logP`: bit-length of the special modulus.
- `logQ`: bit-length of the ciphertext modulus.
- `dlen`: how many prime moduli we want to gather for the gadget decomposition. If 0, it will be calculated automatically.
- `scaling_factor`: scaling factor.
"""
struct CKKSParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64
    scaling_factor::UInt128

    CKKSParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, scaling_factor::Integer; dlen::Int64=0) =
        new(ring_param, logP, logQ, dlen, scaling_factor)
end

"""
    CKKSParameters(sketch::CKKSParamSketch)
    CKKSParameters(ring_param::RingParam, P::Union{Missing,Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, scaling_factor::Integer)

CKKSParameters is a struct for the parameters for the CKKS scheme.
"""
struct CKKSParameters
    ring_param::RingParam
    P::Union{Missing,Vector{UInt64}}
    Q::Vector{UInt64}
    dlen::Int64
    scaling_factor::UInt128

    CKKSParameters(ring_param::RingParam, P::Union{Missing,Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, scaling_factor::Integer) =
        new(ring_param, P, Q, dlen, scaling_factor)

    function CKKSParameters(sketch::CKKSParamSketch)
        ring_param, logP, logQ, dlen, Δ = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen, sketch.scaling_factor

        if logP == 0
            @error "logP must be greater than 0"
        else
            maxbits = min(62, logP - 0.5log2(ring_param.m) - 2) # Minimise the key-switching error.
            Plen = ceil(Int64, logP / 62)
            Qlen = ceil(Int64, logQ / maxbits)

            Qbits = round(Int64, logQ / Qlen)
            Q0bit = logQ - Qbits * (Qlen - 1)

            while Q0bit > maxbits
                Q0bit -= Qbits
                Qlen += 1
            end

            Qprimes = find_prime(ring_param, Qbits, Plen + Qlen + 1)
            Q0 = find_prime(ring_param, Q0bit)[1]

            Pbits = logP / Plen
            Pprimes = find_prime(ring_param, Pbits, Plen + 1)
            filter!(x -> x ≠ Q0, Pprimes)
            P = Pprimes[1:Plen]

            filter!(x -> x ∉ P && x ≠ Q0, Qprimes)
            Q = vcat(Q0, Qprimes[1:Qlen-1])

            if dlen == 0
                dlen = floor(Int64, logP / max(Q0bit, Qbits))
            end
        end

        new(ring_param, P, Q, dlen, Δ)
    end
end

export CKKSParamSketch, CKKSParameters