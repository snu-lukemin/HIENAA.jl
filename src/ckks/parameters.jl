"""
    CKKSParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, scaling_factor::Real; dlen::Int64=0)

CKKSParamSketch is a struct for the parameters for the CKKS scheme.
The following parameters stand for:
- `ring_param`: base ring parameters. Currently, only power-of-two cyclotomic rings are supported.
- `logP`: bit-length of the special modulus.
- `logQ`: bit-length of the ciphertext modulus.
- `dlen`: how many prime moduli we want to gather for the gadget decomposition. If 0, it will be calculated automatically.
- `scaling_factor`: scaling factor.
"""
struct CKKSParamSketch <: HEParamSketch
    ring_param::RingParam
    logP::Int64
    logQ::Int64
    dlen::Int64
    scaling_factor::BigFloat

    CKKSParamSketch(ring_param::RingParam, logP::Int64, logQ::Int64, scaling_factor::Real; dlen::Int64=0)::CKKSParamSketch =
        new(ring_param, logP, logQ, dlen, BigFloat(scaling_factor, precision=192))
end

"""
    CKKSParameters(sketch::CKKSParamSketch)
    CKKSParameters(ring_param::RingParam, P::Union{Missing,Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, scaling_factor::Real)

CKKSParameters is a struct for the parameters for the CKKS scheme.
"""
struct CKKSParameters <: HEParameters
    ring_param::RingParam
    P::Union{Missing,Vector{UInt64}}
    Q::Vector{UInt64}
    dlen::Int64
    scaling_factor::BigFloat

    CKKSParameters(ring_param::RingParam, P::Union{Missing,Vector{UInt64}}, Q::Vector{UInt64}, dlen::Int64, scaling_factor::Real)::CKKSParameters =
        new(ring_param, P, Q, dlen, BigFloat(scaling_factor, precision=192))

    function CKKSParameters(sketch::CKKSParamSketch)::CKKSParameters
        ring_param, logP, logQ, dlen, Δ = sketch.ring_param, sketch.logP, sketch.logQ, sketch.dlen, sketch.scaling_factor

        if logP == 0
            error("logP must be greater than 0")
        else
            maxbits = min(62, logP)
            Plen = ceil(Int64, logP / 62)
            Qlen = ceil(Int64, logQ / maxbits)
            Qbits = logQ / Qlen

            Pbits = logP / Plen
            P = find_prime(ring_param, Pbits, Plen)

            Qprimes = find_prime(ring_param, Qbits, Plen + Qlen)
            filter!(x -> x ∉ P, Qprimes)
            Q = Qprimes[1:Qlen]

            if dlen == 0
                dlen = floor(Int64, logP / Qbits)
            end
        end

        new(ring_param, P, Q, dlen, Δ)
    end
end

struct CKKSBootParameters <: HEBootParameters
    gap::Float64
    packlen::Int64

    P0_bits::Float64
    Q0_bits::Float64
    sparse_hw::Int64

    c2s_radix::Int64
    c2s_bits::Float64

    mod_bits::Float64
    K::Int64
    sine_fold::Int64
    sine_degree::Int64
    inverse_degree::Int64

    s2c_radix::Int64
    s2c_bits::Float64

    isthin::Bool

    function CKKSBootParameters(; gap::Real=32.0, packlen::Int64=0, P0_bits::Real=61.0, Q0_bits::Real=60.0, sparse_hw::Int64=32,
        c2s_radix::Int64=3, c2s_bits::Real=58.0, mod_bits::Real=50.0, K::Int64=16, sine_fold::Int64=3, sine_degree::Int64=57, inverse_degree::Int64=15,
        s2c_radix::Int64=3, s2c_bits::Real=42.0, isthin::Bool=false)::CKKSBootParameters

        if Q0_bits > 62
            throw(DomainError("Q0_bits should be less than 62."))
        end
        
        new(Float64(gap), packlen, Float64(P0_bits), Float64(Q0_bits), sparse_hw, c2s_radix, Float64(c2s_bits), Float64(mod_bits), K, sine_fold, sine_degree, inverse_degree, 
        s2c_radix, Float64(s2c_bits), isthin)
    end
end