"""
SubringNTTransformer is a struct which support the number theoretic transform (NTT) over the subring.
You can generate a SubringNTTransformer over a ring Z[X]/(Φₘ(X), Q) with length N = (m-1)/d by ```SubringNTTransformer(m, d, Q)```
"""
struct SubringNTTransformer <: NTTransformer
    Q::Modulus
    m::Int64
    d::UInt64
    N::Int64
    g::Int64
    minv::UInt64
    Ψ::Vector{UInt64}
    Ψinv::Vector{UInt64}
    buff::Vector{UInt64}
    gpowN::Vector{Int32}
    ginvpowN::Vector{Int32}
    ntter::CyclicNTTransformer2a3b5c7d

    function SubringNTTransformer(m::Int64, d::Int64, Q::Modulus)::SubringNTTransformer
        if !isprime(m)
            throw(DomainError("$(m) is not a prime number."))
        end
        if (m - 1) % d ≠ 0
            throw(DomainError("$(d) should divide $(m-1)."))
        end
        if Q.Q % m ≠ 1
            throw(DomainError("$(Q.Q) does not have enough factors to perform NTT."))
        end

        g = primitive_root_finder(m)
        N = (m - 1) ÷ d

        convlen = is2a3b5c7d(N) ? N : next2a3b5c7d(2N - 1)
        ntter = CyclicNTTransformer2a3b5c7d(convlen, Q)

        ξ = primitive_root_finder(Q)
        ζ = ith_root_from_primitive_root(m, ξ, Q)
        ζinv = powermod(ζ, m - 1, Q)

        Ψ, Ψinv = zeros(UInt64, convlen), zeros(UInt64, convlen)
        gᴺ = powermod(g, N, m)
        ζgʲᴺ, ζinvgʲᴺ = ζ, ζinv
        @inbounds for _ = 1:d
            ζgʲᴺ, ζinvgʲᴺ = powermod(ζgʲᴺ, gᴺ, Q), powermod(ζinvgʲᴺ, gᴺ, Q)
            ζgʲᴺgⁱ, ζinvgʲᴺgⁱ = ζgʲᴺ, ζinvgʲᴺ
            for i = 1:N
                Ψ[i] = add(Ψ[i], ζgʲᴺgⁱ, Q)
                Ψinv[i] = add(Ψinv[i], ζinvgʲᴺgⁱ, Q)
                ζgʲᴺgⁱ, ζinvgʲᴺgⁱ = powermod(ζgʲᴺgⁱ, g, Q), powermod(ζinvgʲᴺgⁱ, g, Q)
            end
        end

        ntt!(Ψ, ntter)
        ntt!(Ψinv, ntter)

        buff = Vector{UInt64}(undef, convlen)

        ginv = invmod(g, m)
        gpowN = Int32[powermod(g, i, m) for i = 0:N-1]
        ginvpowN = Int32[powermod(ginv, i, m) for i = 0:N-1]

        new(Q, m, d, N, g, invmod(m, Q), Ψ, Ψinv, buff, gpowN, ginvpowN, ntter)
    end
end

@views function ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::SubringNTTransformer)::Nothing
    buff, N, Q = ntter.buff, ntter.N, ntter.Q

    buff[1] = a[1]
    @. buff[2:N] = a[N:-1:2]
    @. buff[N+1:end] = 0

    ntt!(buff, ntter.ntter)
    lazy_Bmul_to!(buff, buff, ntter.Ψ, Q)
    intt!(buff, ntter.ntter)

    length(buff) ≠ N && add_to!(buff[1:N-1], buff[1:N-1], buff[N+1:2N-1], Q)
    copy!(res, buff[1:N])

    return nothing
end

ntt!(a::AbstractVector{UInt64}, ntter::SubringNTTransformer)::Nothing = ntt_to!(a, a, ntter)

@views function intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::SubringNTTransformer)::Nothing
    buff, N, d, minv, Q = ntter.buff, ntter.N, ntter.d, ntter.minv, ntter.Q

    buff[1] = a[1]
    @. buff[2:N] = a[N:-1:2]
    @. buff[N+1:end] = 0

    ntt!(buff, ntter.ntter)
    lazy_Bmul_to!(buff, buff, ntter.Ψinv, ntter.Q)
    intt!(buff, ntter.ntter)

    length(buff) ≠ N && add_to!(buff[1:N-1], buff[1:N-1], buff[N+1:2N-1], Q)

    t = UInt64(0)
    @inbounds for i = 1:N
        t = add(t, a[i], Q)
    end
    sub_to!(buff[1:N], buff[1:N], Bmul(t, d, Q), Q)
    Bmul_to!(buff[1:N], minv, buff[1:N], Q)
    copy!(res, buff[1:N])

    return nothing
end

intt!(a::AbstractVector{UInt64}, ntter::SubringNTTransformer)::Nothing = intt_to!(a, a, ntter)