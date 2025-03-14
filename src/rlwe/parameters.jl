"""
RLWEParamSketch is a struct for a very vague parameter setting.
"""
struct RLWEParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64

    RLWEParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, dlen::Int64=0)::RLWEParamSketch = 
        new(ring_param, logP, logQ, dlen)
end

"""
RLWEParameters is a struct for a concrete parameter setting.
"""
struct RLWEParameters
    ring_param::RingParam
    P::Union{Missing,Vector{UInt64}}
    Q::Vector{UInt64}
    dlen::Int64

    RLWEParameters(ring_param::RingParam, P::Union{Missing, Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64)::RLWEParameters = 
        new(ring_param, P, Q, dlen)

    function RLWEParameters(sketch::RLWEParamSketch)::RLWEParameters
        ring_param, logP, logQ, dlen = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen
    
        if logP == 0
            P = missing
    
            Qlen = ceil(Int64, logQ / 62)
            Qbits = logQ / Qlen
            Qprimes = find_prime(ring_param, Qbits, Qlen)
            Q = Qprimes
    
            if dlen == 0
                dlen = 1
            end
        else
            maxbits = min(62, logP - 0.5log2(ring_param.m) - 2) # Minimise the key-switching error.
            Plen = ceil(Int64, logP / 62)
            Qlen = ceil(Int64, logQ / maxbits)

            Qbits = logQ / Qlen
            Qprimes = find_prime(ring_param, Qbits, Plen + Qlen + 1)
            
            Pbits = logP / Plen
            P = find_prime(ring_param, Pbits, Plen)

            filter!(x -> x ∉ P, Qprimes)
            Q = Qprimes[1:Qlen]

            if dlen == 0
                dlen = floor(Int64, logP / Qbits)
            end
        end
    
        new(ring_param, P, Q, dlen)
    end
end