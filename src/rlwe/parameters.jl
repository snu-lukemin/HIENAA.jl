"""
RLWEParamSketch is a struct for a very vague parameter setting.
"""
struct RLWEParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64

    RLWEParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, dlen::Int64=0) = new(ring_param, logP, logQ, dlen)
end

"""
RLWEParameters is a struct for a concrete parameter setting.
"""
struct RLWEParameters
    ring_param::RingParam
    P::Union{Missing,Vector{UInt64}}
    Q::Vector{UInt64}
    dlen::Int64

    RLWEParameters(ring_param::RingParam, P::Union{Missing, Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64) = 
        new(ring_param, P, Q, dlen)

    function RLWEParameters(sketch::RLWEParamSketch)
        ring_param, logP, logQ, dlen = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen
    
        if logP == 0
            P = missing
    
            Qlen = ceil(Int64, logQ / 62)
            Qbits = round(Int64, logQ / Qlen)
            Q0bit = logQ - Qbits * (Qlen - 1)
    
            while Q0bit > 62
                Q0bit -= Qbits
                Qlen += 1
            end
    
            Qprimes = find_prime(ring_param, Qbits, Qlen)
            Q0 = find_prime(ring_param, Q0bit)[1]
            Q = Q0bit == Qbits ? Qprimes : vcat(Q0, Qprimes[1:end-1])
    
            if dlen == 0
                dlen = 1
            end
        else
            Plen = ceil(Int64, logP / 62)
            Pbits = logP / Plen
            P = find_prime(ring_param, Pbits, Plen)
    
            Qlen = ceil(Int64, logQ / 62)
            Qbits = round(Int64, logQ / Qlen)
            Q0bit = logQ - Qbits * (Qlen - 1)
    
            while Q0bit > 62
                Q0bit -= Qbits
                Qlen += 1
            end
    
            if isapprox(Pbits, Qbits, atol=0.005)
                Qprimes = find_prime(ring_param, Qbits, Plen + Qlen)
                Q0 = find_prime(ring_param, Q0bit)[1]
                Q = (Q0bit == Qbits) ? Qprimes[Plen+1:end] : vcat(Q0, Qprimes[Plen+1:end-1])
            else
                Qprimes = find_prime(ring_param, Qbits, Qlen)
                Q0 = find_prime(ring_param, Q0bit)[1]
                Q = (Q0bit == Qbits) ? Qprimes : vcat(Q0, Qprimes[1:end-1])
            end
    
            if dlen == 0
                dlen = floor(Int64, logP / max(Q0bit, Qbits))
            end
        end
    
        new(ring_param, P, Q, dlen)
    end
end

export RLWEParamSketch, RLWEParameters