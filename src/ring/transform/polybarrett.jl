abstract type Reductor end
abstract type ReductorCyclo <: Reductor end
abstract type ReductorArb <: Reductor end

struct ReductorArbNTT <: ReductorArb
    deg::Int64
    N::Int64
    Q::Modulus
    α::Int64
    A::Int64
    ñ::Int64
    ntterñ::CyclicNTTransformer2a3b5c7d
    ntterA::CyclicNTTransformer2a3b5c7d
    poly::Vector{UInt64}
    Xdeg_over_poly::Vector{UInt64}
    fbuff::Vector{UInt64}
    rbuff::Vector{UInt64}

    @views function ReductorArbNTT(deg::Int64, poly::Vector{Int64}, Q::Modulus)::ReductorArbNTT
        N = length(poly) - 1

        ñ = next2a3b(N)

        α = deg - N
        A = next2a3b(2α + 1)

        ntterñ = CyclicNTTransformer2a3b5c7d(ñ, Q)
        ntterA = CyclicNTTransformer2a3b5c7d(A, Q)

        poly = Bred.(poly, Ref(Q))
        Xdeg_over_poly = division_mod_Q(vcat(zeros(UInt64, deg), UInt64(1)), poly, Q)

        poly = zeropadto(poly, ñ)
        Xdeg_over_poly = zeropadto(Xdeg_over_poly, A)

        ntt!(poly, ntterñ)
        ntt!(Xdeg_over_poly, ntterA)

        fbuff = Vector{UInt64}(undef, A)
        rbuff = Vector{UInt64}(undef, ñ)

        new(deg, N, Q, α, A, ñ, ntterñ, ntterA, poly, Xdeg_over_poly, fbuff, rbuff)
    end
end

# We use Optimised Barrett reduction for polynomial, from https://eprint.iacr.org/2017/748.
# In this implementation, Qₛₚ = (Xᵐ-1)/(X^(m/p)-1) for the smallest prime factor p.
# In the original paper, the authors used polynomials with smaller degree than ours.
# However, for a generalised code, we shall stick to these polynomials.
struct ReductorCycloNTT <: ReductorCyclo
    p::Int64
    moverp::Int64
    iseasy::Bool
    m::Int64
    deg::Int64
    N::Int64
    Q::Modulus
    α::Int64
    A::Int64
    ñ::Int64
    ntterñ::Union{CyclicNTTransformer2a3b5c7d,Missing}
    ntterA::Union{CyclicNTTransformer2a3b5c7d,Missing}
    poly::Vector{UInt64}
    Xdeg_over_poly::Vector{UInt64}
    fbuff::Vector{UInt64}
    rbuff::Vector{UInt64}

    @views function ReductorCycloNTT(m::Int64, Q::Modulus)::ReductorCycloNTT
        N = totient(m)
        factors = factor(Vector, m)
        p = factors[1]
        moverp = m ÷ p
        deg = m - moverp

        if deg == N
            new(p, moverp, true, m, deg, N, Q, 0, 0, 0, missing, missing, UInt64[], UInt64[], UInt64[], UInt64[])
        else
            ñ = next2a3b(N)

            α = deg - N
            A = next2a3b(2α + 1)

            ntterñ = CyclicNTTransformer2a3b5c7d(ñ, Q)
            ntterA = CyclicNTTransformer2a3b5c7d(A, Q)

            poly = cyclotomic_finder(m)
            Xdeg_over_poly = division(vcat(zeros(Int64, deg), 1), poly)

            poly = Bred.(zeropadto(poly, ñ), Ref(Q))
            Xdeg_over_poly = Bred.(zeropadto(Xdeg_over_poly, A), Ref(Q))

            ntt!(poly, ntterñ)
            ntt!(Xdeg_over_poly, ntterA)

            fbuff = Vector{UInt64}(undef, A)
            rbuff = Vector{UInt64}(undef, ñ)

            new(p, moverp, false, m, deg, N, Q, α, A, ñ, ntterñ, ntterA, poly, Xdeg_over_poly, fbuff, rbuff)
        end
    end
end

"""
barrett! computes a (mod poly).
"""
@views function barrett!(a::AbstractVector{UInt64}, rdtor::Union{ReductorArbNTT,ReductorCycloNTT})
    deg, N, α, ñ, Q = rdtor.deg, rdtor.N, rdtor.α, rdtor.ñ, rdtor.Q

    # Compute f = ⌊a/Xᴺ⌋
    @inbounds @simd for i = 1:α
        rdtor.fbuff[i] = a[N+i]
    end
    @. rdtor.fbuff[α+1:end] = zero(UInt64)

    # Compute f = f × ⌊Xᴺ⁺ᵅ/p⌋
    ntt!(rdtor.fbuff, rdtor.ntterA)
    lazy_Bmul_to!(rdtor.fbuff, rdtor.fbuff, rdtor.Xdeg_over_poly, Q)
    intt!(rdtor.fbuff, rdtor.ntterA)

    # Compute r = ⌊f/Xᵅ⌋ % (Xⁿ̃ - 1)
    @inbounds for i = 1:ceil(Int64, α / ñ), j = 1:ñ
        i * ñ + j > α && break
        rdtor.fbuff[α+j] = add(rdtor.fbuff[α+j], rdtor.fbuff[α+i*ñ+j], Q)
        rdtor.buff[α+i*ñ+j] = 0
    end

    idx = min(α, ñ)
    @inbounds @simd for i = 1:idx
        rdtor.rbuff[i] = rdtor.fbuff[α+i]
    end
    @. rdtor.rbuff[idx+1:end] = zero(UInt64)

    # Compute r = r × p (mod Xⁿ̃ - 1)
    ntt!(rdtor.rbuff, rdtor.ntterñ)
    lazy_Bmul_to!(rdtor.rbuff, rdtor.rbuff, rdtor.poly, Q)
    intt!(rdtor.rbuff, rdtor.ntterñ)

    # Compute a = a % (Xⁿ̃ - 1)
    @inbounds for i = 1:ceil(Int64, deg / ñ), j = 1:ñ
        i * ñ + j > deg && break
        a[j] = add(a[j], a[i*ñ+j], Q)
        a[i*ñ+j] = 0
    end

    # Compute a -= r
    @inbounds for i = 1:N
        a[i] = sub(a[i], rdtor.rbuff[i], Q)
    end
    @inbounds @simd for i = N+1:length(a)
        a[i] = zero(UInt64)
    end

    return nothing
end

#===============================================================================#

"""
ReductorArbWord
"""
struct ReductorArbWord <: ReductorArb
    deg::Int64
    N::Int64
    Q::Modulus
    Plen::Int64
    P::Vector{Modulus}
    α::Int64
    A::Int64
    ñ::Int64
    beP2Q::BasisExtender
    ntterñ::Vector{CyclicNTTransformer2a3b5c7d}
    ntterA::Vector{CyclicNTTransformer2a3b5c7d}
    poly::Vector{Vector{UInt64}}
    Xdeg_over_poly::Vector{Vector{UInt64}}
    fbuff::Vector{Vector{UInt64}}
    rbuff::Vector{Vector{UInt64}}
    buff::Vector{Vector{UInt64}}

    @views function ReductorArbWord(deg::Int64, poly::Vector{Int64}, Q::Modulus)
        N = length(poly) - 1

        ñ = next2a3b(N)
        α = deg - N
        A = next2a3b(2α + 1)

        Plen = ceil(Int64, (2log2(Q.Q) + log2(max(ñ, A))) / 62)
        P = Modulus.(collect(Ring.find_prime_reductor(deg, N, 62, Plen)))

        beP2Q = BasisExtender(P, [Q])

        ntterñ = CyclicNTTransformer2a3b5c7d[CyclicNTTransformer2a3b5c7d(ñ, Pi) for Pi = P]
        ntterA = CyclicNTTransformer2a3b5c7d[CyclicNTTransformer2a3b5c7d(A, Pi) for Pi = P]

        _poly = Bred.(poly, Ref(Q))
        _Xdeg_over_poly = division_mod_Q(vcat(zeros(UInt64, N + α), UInt64(1)), _poly, Q)

        poly = Vector{Vector{UInt64}}(undef, Plen)
        Xdeg_over_poly = Vector{Vector{UInt64}}(undef, Plen)

        @inbounds for i = 1:Plen
            poly[i] = zeros(UInt64, ñ)
            Xdeg_over_poly[i] = zeros(UInt64, A)

            for j = eachindex(_poly)
                poly[i][j] = Bred(_poly[j], P[i])
            end
            for j = eachindex(_Xdeg_over_poly)
                Xdeg_over_poly[i][j] = Bred(_Xdeg_over_poly[j], P[i])
            end

            ntt!(poly[i], ntterñ[i])
            ntt!(Xdeg_over_poly[i], ntterA[i])
        end

        fbuff = [Vector{UInt64}(undef, A) for _ = 1:Plen]
        rbuff = [Vector{UInt64}(undef, ñ) for _ = 1:Plen]
        buff = [Vector{UInt64}(undef, max(A, min(deg, ñ)))]

        new(deg, N, Q, Plen, P, α, A, ñ, beP2Q, ntterñ, ntterA, poly, Xdeg_over_poly, fbuff, rbuff, buff)
    end
end

struct ReductorCycloWord <: ReductorCyclo
    p::Int64
    moverp::Int64
    iseasy::Bool
    m::Int64
    deg::Int64
    N::Int64
    Q::Modulus
    Plen::Int64
    P::Vector{Modulus}
    α::Int64
    A::Int64
    ñ::Int64
    beP2Q::Union{BasisExtender,Missing}
    ntterñ::Vector{CyclicNTTransformer2a3b5c7d}
    ntterA::Vector{CyclicNTTransformer2a3b5c7d}
    poly::Vector{Vector{UInt64}}
    Xdeg_over_poly::Vector{Vector{UInt64}}
    fbuff::Vector{Vector{UInt64}}
    rbuff::Vector{Vector{UInt64}}
    buff::Vector{Vector{UInt64}}

    function ReductorCycloWord(m::Int64, Q::Modulus)::ReductorCycloWord
        N = totient(m)
        factors = factor(Vector, m)
        p = factors[1]
        moverp = m ÷ p
        deg = m - moverp

        if deg == N
            new(p, moverp, true, m, deg, N, Q, 0, Modulus[], 0, 0, 0, missing, CyclicNTTransformer2a3b5c7d[], CyclicNTTransformer2a3b5c7d[], Vector{UInt64}[], Vector{UInt64}[], Vector{UInt64}[], Vector{UInt64}[], Vector{UInt64}[])
        else
            ñ = next2a3b(N)
            α = deg - N
            A = next2a3b(2α + 1)

            Plen = ceil(Int64, (2log2(Q.Q) + log2(max(ñ, A))) / 62)
            P = Modulus.(collect(Ring.find_prime_reductor(deg, N, 62, 3)))

            beP2Q = BasisExtender(P[1:Plen], [Q])

            ntterñ = CyclicNTTransformer2a3b5c7d[CyclicNTTransformer2a3b5c7d(ñ, Pi) for Pi = P]
            ntterA = CyclicNTTransformer2a3b5c7d[CyclicNTTransformer2a3b5c7d(A, Pi) for Pi = P]

            _poly = cyclotomic_finder(m)
            _Xdeg_over_poly = division(vcat(zeros(Int64, N + α), 1), _poly)

            poly = Vector{Vector{UInt64}}(undef, 3)
            Xdeg_over_poly = Vector{Vector{UInt64}}(undef, 3)

            @inbounds for i = 1:3
                poly[i] = zeros(UInt64, ñ)
                Xdeg_over_poly[i] = zeros(UInt64, A)

                for j = eachindex(_poly)
                    poly[i][j] = Bred(_poly[j], P[i])
                end
                for j = eachindex(_Xdeg_over_poly)
                    Xdeg_over_poly[i][j] = Bred(_Xdeg_over_poly[j], P[i])
                end

                ntt!(poly[i], ntterñ[i])
                ntt!(Xdeg_over_poly[i], ntterA[i])
            end

            fbuff = [Vector{UInt64}(undef, A) for _ = 1:3]
            rbuff = [Vector{UInt64}(undef, ñ) for _ = 1:3]
            buff = [Vector{UInt64}(undef, max(A, min(deg, ñ)))]

            new(p, moverp, false, m, deg, N, Q, Plen, P, α, A, ñ, beP2Q, ntterñ, ntterA, poly, Xdeg_over_poly, fbuff, rbuff, buff)
        end
    end

    # Update the modulus Q.
    function ReductorCycloWord(rdtor::ReductorCycloWord, Q::Modulus)::ReductorCycloWord
        p, moverp, iseasy, m, deg, N, Plen, P, α, A, ñ, beP2Q, ntterñ, ntterA, poly, Xdeg_over_poly, fbuff, rbuff, buff =
            rdtor.p, rdtor.moverp, rdtor.iseasy, rdtor.m, rdtor.deg, rdtor.N, rdtor.Plen, rdtor.P, rdtor.α, rdtor.A, rdtor.ñ, rdtor.beP2Q, rdtor.ntterñ, rdtor.ntterA, rdtor.poly, rdtor.Xdeg_over_poly, rdtor.fbuff, rdtor.rbuff, rdtor.buff

        if !iseasy
            Plen = ceil(Int64, (2log2(Q.Q) + log2(max(ñ, A))) / 62)
            beP2Q = BasisExtender(P[1:Plen], [Q])
        end

        new(p, moverp, iseasy, m, deg, N, Q, Plen, P, α, A, ñ, beP2Q, ntterñ, ntterA, poly, Xdeg_over_poly, fbuff, rbuff, buff)
    end
end

@views function barrett!(a::AbstractVector{UInt64}, rdtor::Union{ReductorArbWord,ReductorCycloWord})::Nothing
    deg, N, α, ñ, Q, Plen, P, beP2Q = rdtor.deg, rdtor.N, rdtor.α, rdtor.ñ, rdtor.Q, rdtor.Plen, rdtor.P, rdtor.beP2Q

    # Compute f = ⌊a/Xᴺ⌋
    @inbounds for j = 1:Plen
        @. rdtor.fbuff[j][1:α] = a[N+1:N+α]
        @. rdtor.fbuff[j][α+1:end] = zero(UInt64)
    end

    # Compute f = f × ⌊Xᴺ⁺ᵅ/p⌋
    @inbounds for j = 1:Plen
        ntt!(rdtor.fbuff[j], rdtor.ntterA[j])
        lazy_Bmul_to!(rdtor.fbuff[j], rdtor.fbuff[j], rdtor.Xdeg_over_poly[j], P[j])
        intt!(rdtor.fbuff[j], rdtor.ntterA[j])
    end
    basis_extend_to!(rdtor.buff, 1:α, rdtor.fbuff, α+1:2α, beP2Q)

    # Compute r = ⌊f/Xᵅ⌋ % (Xⁿ̃ - 1)
    @inbounds for i = 1:ceil(Int64, α / ñ), j = 1:ñ
        i * ñ + j > α && break
        rdtor.buff[1][j] = add(rdtor.buff[1][j], rdtor.buff[1][i*ñ+j], Q)
        rdtor.buff[1][i*ñ+j] = 0
    end

    idx = min(α, ñ)
    @inbounds for j = 1:Plen
        @. rdtor.rbuff[j][1:idx] = rdtor.buff[1][1:idx]
        @. rdtor.rbuff[j][idx+1:end] = zero(UInt64)
    end

    # Compute r = r × p (mod Xⁿ̃ - 1)
    idx = min(deg, ñ)
    @inbounds for j = 1:Plen
        ntt!(rdtor.rbuff[j], rdtor.ntterñ[j])
        lazy_Bmul_to!(rdtor.rbuff[j], rdtor.rbuff[j], rdtor.poly[j], P[j])
        intt!(rdtor.rbuff[j], rdtor.ntterñ[j])
    end
    basis_extend_to!(rdtor.buff, 1:idx, rdtor.rbuff, 1:idx, beP2Q)

    # Compute a = a % (Xⁿ̃ - 1)
    @inbounds for i = 1:ceil(Int64, deg / ñ), j = 1:ñ
        i * ñ + j > deg && break
        a[j] = add(a[j], a[i*ñ+j], Q)
        a[i*ñ+j] = 0
    end

    # Compute a -= r
    @inbounds for i = 1:N
        a[i] = sub(a[i], rdtor.buff[1][i], Q)
    end
    @inbounds @simd for i = N+1:length(a)
        a[i] = 0
    end

    return nothing
end

"""
reduce! computes a (mod Φₘ).
"""
function reduce!(a::AbstractVector{UInt64}, rdtor::ReductorCyclo)::Nothing
    p, moverp, m, Q = rdtor.p, rdtor.moverp, rdtor.m, rdtor.Q

    # reduction by Qₛₚ
    @inbounds for j = 1:moverp
        for i = 0:p-2
            a[i*moverp+j] = sub(a[i*moverp+j], a[m-moverp+j], Q)
        end
        a[m-moverp+j] = 0
    end

    !rdtor.iseasy && barrett!(a, rdtor)

    return nothing
end

reduce!(a::AbstractVector{UInt64}, rdtor::ReductorArb)::Nothing = barrett!(a, rdtor)