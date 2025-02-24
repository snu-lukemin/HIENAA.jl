"""
    BGVParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, ptxt_modulus::Int64; dlen::Int64=0, ispacking::Bool=true)

BGVParamSketch is a struct for a very vague parameter setting for the BGV scheme.
The following parameters stand for:
- `ring_param`: base ring parameters.
- `logP`: bit-length of the special modulus.
- `logQ`: bit-length of the ciphertext modulus.
- `dlen`: how many prime moduli we want to gather for the gadget decomposition. If 0, it will be calculated automatically.
- `ptxt_modulus`: plaintext modulus.
- `ispacking`: whether the messages are packed or not.
"""
struct BGVParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64
    ptxt_modulus::Int64
    ispacking::Bool

    BGVParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, ptxt_modulus::Int64; dlen::Int64=0, ispacking::Bool=true) =
        new(ring_param, logP, logQ, dlen, ptxt_modulus, ispacking)
end

"""
    BGVParameters(sketch::BGVParamSketch)
    BGVParameters(ring_param::RingParam, P::Union{Missing,Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, ptxt_modulus::Int64, ispacking::Bool)

BGVParameters is a struct for the parameters for the BGV scheme.
"""
struct BGVParameters
    ring_param::RingParam
    P::Union{Missing,Vector{UInt64}}
    Q::Vector{UInt64}
    dlen::Int64
    ptxt_modulus::UInt64
    ispacking::Bool

    BGVParameters(ring_param::RingParam, P::Union{Missing,Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, ptxt_modulus::Int64, ispacking::Bool) =
        new(ring_param, P, Q, dlen, ptxt_modulus, ispacking)

    function BGVParameters(sketch::BGVParamSketch)
        ring_param, logP, logQ, dlen, t, ispacking = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen, sketch.ptxt_modulus, sketch.ispacking

        if logP == 0
            @error "logP must be greater than 0"
        else
            maxbits = min(62, logP)
            Plen = ceil(Int64, logP / 62)
            Qlen = ceil(Int64, logQ / maxbits)
            
            Qbits = round(Int64, logQ / Qlen)
            Q0bit = logQ - Qbits * (Qlen - 1)

            while Q0bit > maxbits
                Q0bit -= Qbits
                Qlen += 1
            end

            Qprimes = find_prime(ring_param, Qbits, Plen + Qlen + 1, onemod=t)
            Q0 = find_prime(ring_param, Q0bit, onemod=t)[1]

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

        new(ring_param, P, Q, dlen, t, ispacking)
    end
end

export BGVParamSketch, BGVParameters