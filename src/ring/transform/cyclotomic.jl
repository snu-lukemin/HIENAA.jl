abstract type CyclotomicNTTransformer <: NTTransformer end

struct CyclotomicNTTransformerPow2 <: CyclotomicNTTransformer
    Q::Modulus
    m::Int64
    N::Int64
    N⁻¹::UInt64
    Ψ::Vector{UInt64}
    Ψinv::Vector{UInt64}

    function CyclotomicNTTransformerPow2(m::Int64, Q::Modulus)::CyclotomicNTTransformerPow2
        if Q.Q & (m - 1) ≠ 1
            throw(DomainError("ϕ(Q) does not have enough 2 factor to perform NTT."))
        end

        N = m >> 1

        ξ = primitive_root_finder(Q)
        ζ = ith_root_from_primitive_root(m, ξ, Q)
        ζinv = powermod(ζ, 2N - 1, Q)
        Ψ, Ψinv = Vector{UInt64}(undef, N), Vector{UInt64}(undef, N)
        Ψ[1], Ψinv[1], Ψ[2], Ψinv[2] = UInt64(1), UInt64(1), ζ, ζinv
        @inbounds for i = 3:N
            Ψ[i], Ψinv[i] = Bmul(ζ, Ψ[i-1], Q), Bmul(ζinv, Ψinv[i-1], Q)
        end
        scramble!(Ψ, 2)
        scramble!(Ψinv, 2)

        Mform!(Ψ, Q)
        Mform!(Ψinv, Q)

        new(Q, m, N, invmod(N, Q), Ψ, Ψinv)
    end
end

@views ntt!(a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerPow2)::Nothing = begin
    _ntt_2a3b5c7d!(a, ntter.Ψ, ntter.Ψ, ntter.Ψ, ntter.Ψ, ntter.Q)

    return nothing
end

@views ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerPow2)::Nothing = begin
    copy!(res, a)
    ntt!(res, ntter)

    return nothing
end

@views intt!(a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerPow2)::Nothing = begin
    _intt_2a3b5c7d!(a, ntter.Ψinv, ntter.Ψinv, ntter.Ψinv, ntter.Ψinv, ntter.Q)
    Bmul_to!(a, ntter.N⁻¹, a, ntter.Q)

    return nothing
end

@views intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerPow2)::Nothing = begin
    copy!(res, a)
    intt!(res, ntter)

    return nothing
end

# We use Bluestein NTT for the sake of the low space complexity. 
# This approach is somewhat similar to HElib, although we simply apply Bluestein NTT of degree m.
# In the original cuHE and polynomial Barrett reduction paper, the authors make use of power-of-two NTT instead.
# Although the time complexity is at least half, the space complexity of their approach is as twice as ours.
struct CyclotomicNTTransformerArb <: CyclotomicNTTransformer
    Q::Modulus
    m::Int64
    N::Int64
    ntter::CyclicNTTransformerBluestein
    rdtor::ReductorCycloNTT
    buff::Vector{UInt64}
    nttidx::Vector{Int64}
    autidxset::Vector{NTuple{Num,Int64}} where {Num}
    dims::Vector{Int64}
    gens::Vector{Int64}

    function CyclotomicNTTransformerArb(m::Int64, Q::Modulus)::CyclotomicNTTransformerArb
        if iseven(m)
            throw(DomainError("Even cyclotomic degree is not supported, unless power of two."))
        end

        dims, gens = find_generators_mod_m(m)
        if length(dims) > 4
            throw(DomainError("Cyclotomic degree with factors more than five is not supported."))
        end

        ntter = CyclicNTTransformerBluestein(m, Q)
        rdtor = ReductorCycloNTT(m, Q)
        buff = Vector{UInt64}(undef, m)

        if length(dims) == 1
            nttidx = Vector{Int64}(undef, dims[1])
            autidxset = Vector{NTuple{1,Int64}}(undef, ntter.N)
            @inbounds for i1 = 0:dims[1]-1
                idx = 1 + i1
                nttidx[idx] = (powermod(gens[1], i1, m)) % m + 1
                autidxset[nttidx[idx]] = (i1,)
            end
        elseif length(dims) == 2
            nttidx = Vector{Int64}(undef, dims[1] * dims[2])
            autidxset = Vector{NTuple{2,Int64}}(undef, ntter.N)
            @inbounds for i1 = 0:dims[1]-1, i2 = 0:dims[2]-1
                idx = 1 + i1 + dims[1] * i2
                nttidx[idx] = (powermod(gens[1], i1, m) * powermod(gens[2], i2, m)) % m + 1
                autidxset[nttidx[idx]] = (i1, i2)
            end
        elseif length(dims) == 3
            nttidx = Vector{Int64}(undef, dims[1] * dims[2] * dims[3])
            autidxset = Vector{NTuple{3,Int64}}(undef, ntter.N)
            @inbounds for i1 = 0:dims[1]-1, i2 = 0:dims[2]-1, i3 = 0:dims[3]-1
                idx = 1 + i1 + dims[1] * i2 + dims[1] * dims[2] * i3
                nttidx[idx] = (powermod(gens[1], i1, m) * powermod(gens[2], i2, m) * powermod(gens[3], i3, m)) % m + 1
                autidxset[nttidx[idx]] = (i1, i2, i3)
            end
        else
            nttidx = Vector{Int64}(undef, dims[1] * dims[2] * dims[3] * dims[4])
            autidxset = Vector{NTuple{4,Int64}}(undef, ntter.N)
            @inbounds for i1 = 0:dims[1]-1, i2 = 0:dims[2]-1, i3 = 0:dims[3]-1, i4 = 0:dims[4]-1
                idx = 1 + i1 + dims[1] * i2 + dims[1] * dims[2] * i3 + dims[1] * dims[2] * dims[3] * i4
                nttidx[idx] = (powermod(gens[1], i1, m) * powermod(gens[2], i2, m) * powermod(gens[3], i3, m) * powermod(gens[4], i4, m)) % m + 1
                autidxset[nttidx[idx]] = (i1, i2, i3, i4)
            end
        end

        new(Q, m, rdtor.N, ntter, rdtor, buff, nttidx, autidxset, dims, gens)
    end
end

# Bluestein NTT
@views function ntt!(a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerArb)::Nothing
    @. ntter.buff = UInt64(0)
    @inbounds @simd for i = 1:ntter.N
        ntter.buff[i] = a[i]
    end

    ntt!(ntter.buff, ntter.ntter)
    @inbounds @simd for i = 1:ntter.N
        a[i] = ntter.buff[ntter.nttidx[i]]
    end

    return nothing
end

@views ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerArb)::Nothing = begin
    copy!(res, a)
    ntt!(res, ntter)

    return nothing
end

# Bluestein iNTT
@views function intt!(a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerArb, islazy::Bool=false)::Nothing
    @. ntter.buff = 0
    @inbounds @simd for i = 1:ntter.N
        ntter.buff[ntter.nttidx[i]] = a[i]
    end

    intt!(ntter.buff, ntter.ntter)

    if islazy
        if length(a) < ntter.m
            throw(DomainError("The vector length does not match."))
        end
        copy!(a, ntter.buff)
    else
        reduce!(ntter.buff, ntter.rdtor)
        @inbounds @simd for i = 1:ntter.N
            a[i] = ntter.buff[i]
        end
    end

    return nothing
end

@views intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclotomicNTTransformerArb, islazy::Bool=false)::Nothing = begin
    @inbounds @simd for i = eachindex(res)
        res[i] = a[i]
    end
    intt!(res, ntter, islazy)

    return nothing
end

(::Type{CyclotomicNTTransformer})(m::Int64, Q::Modulus)::CyclotomicNTTransformer =
    ispow2(m) ? CyclotomicNTTransformerPow2(m, Q) : CyclotomicNTTransformerArb(m, Q)