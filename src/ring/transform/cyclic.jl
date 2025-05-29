abstract type CyclicNTTransformer <: NTTransformer end

function gen_roots(N::Int64, radix::Int64, ξ::UInt64, Q::Modulus)::Tuple{Vector{UInt64}, Vector{UInt64}}
    (N == 1) && return UInt64[1], UInt64[1]

    ζ = ith_root_from_primitive_root(N, ξ, Q)
    ζinv = powermod(ζ, N - 1, Q)

    Ψ, Ψinv = Vector{UInt64}(undef, N ÷ radix), Vector{UInt64}(undef, N ÷ radix)
    Ψ[begin], Ψinv[begin] = UInt64(1), UInt64(1)
    @inbounds for i = 2:N÷radix
        Ψ[i], Ψinv[i] = Bmul(ζ, Ψ[i-1], Q), Bmul(ζinv, Ψinv[i-1], Q)
    end

    scramble!(Ψ, radix)
    scramble!(Ψinv, radix)
    resize!(Ψ, N + radix - 1)
    resize!(Ψinv, N + radix - 1)

    len = N ÷ radix
    @inbounds for i = 1:len
        Ψ[i+len] = Ψ[i]
        Ψinv[i+len] = Ψinv[i]
        for j = 2len:len:(radix-1)*len
            Ψ[i+j] = Bmul(Ψ[i+j-len], Ψ[i], Q)
            Ψinv[i+j] = Bmul(Ψinv[i+j-len], Ψinv[i], Q)
        end
    end

    while len > 1
        len ÷= radix
        @inbounds @simd for i = 1:len
            Ψ[len+i] = Ψ[radix*len+i]
            Ψinv[len+i] = Ψinv[radix*len+i]
        end
        @inbounds for j = 2:radix-1
            for i = 1:len
                Ψ[len*j+i] = Bmul(Ψ[len*(j-1)+i], Ψ[len+i], Q)
                Ψinv[len*j+i] = Bmul(Ψinv[len*(j-1)+i], Ψinv[len+i], Q)
            end
        end
    end

    ω = ith_root_from_primitive_root(radix, ξ, Q)
    Ψ[N+1], Ψinv[N+1] = ω, ω
    @inbounds for i = 2:radix-1
        Ψ[N+i] = Bmul(Ψ[N+i-1], ω, Q)
        Ψinv[N+i] = Ψ[N+i]
    end

    Mform!(Ψ, Q)
    Mform!(Ψinv, Q)

    Ψ, Ψinv
end

"""
An optimised normal input Cooley-Tukey algorithm.
"""
function _ntt_2a3b5c7d!(a::AbstractVector{UInt64}, Ψ2::Vector{UInt64}, Ψ3::Vector{UInt64}, Ψ5::Vector{UInt64}, Ψ7::Vector{UInt64}, Q::Modulus)::Nothing
    N = length(a)
    N2, N3, N5, N7 = factor2357(N)
    N23, N235 = N2 * N3, N2 * N3 * N5

    Mform!(a, Q)
    twoQ = Q.Q << 1
    if N2 > 1
        @inbounds for id = 1:N2:N
            m, logkp1, k = 1, trailing_zeros(N2), N2 >> 1
            for j = id:id+k-1
                t = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                u = lazy_Mmul(a[j+k], Ψ2[2], Q)

                a[j] = t + u
                a[j+k] = t - u + twoQ
            end
            m, logkp1, k = m << 1, logkp1 - 1, k >> 1

            while logkp1 > 3
                for i = 0:m-1
                    j1 = id + i << logkp1
                    j2 = j1 + k - 1
                    for j = j1:8:j2
                        t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                        t1 = a[j+1] ≥ twoQ ? a[j+1] - twoQ : a[j+1]
                        t2 = a[j+2] ≥ twoQ ? a[j+2] - twoQ : a[j+2]
                        t3 = a[j+3] ≥ twoQ ? a[j+3] - twoQ : a[j+3]
                        t4 = a[j+4] ≥ twoQ ? a[j+4] - twoQ : a[j+4]
                        t5 = a[j+5] ≥ twoQ ? a[j+5] - twoQ : a[j+5]
                        t6 = a[j+6] ≥ twoQ ? a[j+6] - twoQ : a[j+6]
                        t7 = a[j+7] ≥ twoQ ? a[j+7] - twoQ : a[j+7]

                        u0 = lazy_Mmul(a[j+k], Ψ2[m+i+1], Q)
                        u1 = lazy_Mmul(a[j+k+1], Ψ2[m+i+1], Q)
                        u2 = lazy_Mmul(a[j+k+2], Ψ2[m+i+1], Q)
                        u3 = lazy_Mmul(a[j+k+3], Ψ2[m+i+1], Q)
                        u4 = lazy_Mmul(a[j+k+4], Ψ2[m+i+1], Q)
                        u5 = lazy_Mmul(a[j+k+5], Ψ2[m+i+1], Q)
                        u6 = lazy_Mmul(a[j+k+6], Ψ2[m+i+1], Q)
                        u7 = lazy_Mmul(a[j+k+7], Ψ2[m+i+1], Q)

                        a[j] = t0 + u0
                        a[j+1] = t1 + u1
                        a[j+2] = t2 + u2
                        a[j+3] = t3 + u3
                        a[j+4] = t4 + u4
                        a[j+5] = t5 + u5
                        a[j+6] = t6 + u6
                        a[j+7] = t7 + u7

                        a[j+k] = t0 - u0 + twoQ
                        a[j+k+1] = t1 - u1 + twoQ
                        a[j+k+2] = t2 - u2 + twoQ
                        a[j+k+3] = t3 - u3 + twoQ
                        a[j+k+4] = t4 - u4 + twoQ
                        a[j+k+5] = t5 - u5 + twoQ
                        a[j+k+6] = t6 - u6 + twoQ
                        a[j+k+7] = t7 - u7 + twoQ
                    end
                end
                m, logkp1, k = m << 1, logkp1 - 1, k >> 1
            end

            if logkp1 == 3
                for i = 0:m-1
                    j = i << 3 + id

                    t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                    t1 = a[j+1] ≥ twoQ ? a[j+1] - twoQ : a[j+1]
                    t2 = a[j+2] ≥ twoQ ? a[j+2] - twoQ : a[j+2]
                    t3 = a[j+3] ≥ twoQ ? a[j+3] - twoQ : a[j+3]

                    u0 = lazy_Mmul(a[j+k], Ψ2[m+i+1], Q)
                    u1 = lazy_Mmul(a[j+k+1], Ψ2[m+i+1], Q)
                    u2 = lazy_Mmul(a[j+k+2], Ψ2[m+i+1], Q)
                    u3 = lazy_Mmul(a[j+k+3], Ψ2[m+i+1], Q)

                    a[j] = t0 + u0
                    a[j+1] = t1 + u1
                    a[j+2] = t2 + u2
                    a[j+3] = t3 + u3

                    a[j+k] = t0 - u0 + twoQ
                    a[j+k+1] = t1 - u1 + twoQ
                    a[j+k+2] = t2 - u2 + twoQ
                    a[j+k+3] = t3 - u3 + twoQ
                end
                m, logkp1, k = m << 1, logkp1 - 1, k >> 1
            end

            if logkp1 == 2
                for i = 0:m-1
                    j = i << 2 + id

                    t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                    t1 = a[j+1] ≥ twoQ ? a[j+1] - twoQ : a[j+1]

                    u0 = lazy_Mmul(a[j+k], Ψ2[m+i+1], Q)
                    u1 = lazy_Mmul(a[j+k+1], Ψ2[m+i+1], Q)

                    a[j] = t0 + u0
                    a[j+1] = t1 + u1

                    a[j+k] = t0 - u0 + twoQ
                    a[j+k+1] = t1 - u1 + twoQ
                end
                m, logkp1, k = m << 1, logkp1 - 1, k >> 1
            end

            if logkp1 == 1
                for i = 0:m-1
                    j = i << 1 + id

                    t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                    u0 = lazy_Mmul(a[j+k], Ψ2[m+i+1], Q)

                    a[j] = t0 + u0
                    a[j+k] = t0 - u0 + twoQ
                end
            end
        end
    end

    if N3 > 1
        ω, ω2 = Ψ3[end-1], Ψ3[end]
        @inbounds for id = 1:N23:N
            m, k = 1, N23 ÷ 3
            while m < N3
                for i = 0:m-1
                    j1 = id + 3k * i
                    j2 = j1 + k - 1
                    ψ1, ψ2 = Ψ3[m+i+1], Ψ3[2m+i+1]
                    for j = j1:j2
                        t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                        t1 = lazy_Mmul(a[j+k], ψ1, Q)
                        t2 = lazy_Mmul(a[j+2k], ψ2, Q)

                        x0 = t0 + t1
                        x0 ≥ twoQ && (x0 -= twoQ)
                        a[j] = x0 + t2
                        t2 = twoQ - t2

                        t0 += t2
                        t1 += t2

                        t0 ≥ twoQ && (t0 -= twoQ)
                        t1 ≥ twoQ && (t1 -= twoQ)

                        x1 = t0 + lazy_Mmul(t1, ω, Q)
                        x2 = t0 + lazy_Mmul(t1, ω2, Q)

                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)

                        a[j+k] = x1
                        a[j+2k] = x2
                    end
                end
                m *= 3
                k ÷= 3
            end
        end
    end

    if N5 > 1
        ω, ω2, ω3, ω4 = Ψ5[end-3], Ψ5[end-2], Ψ5[end-1], Ψ5[end]
        @inbounds for id = 1:N235:N
            m, k = 1, N235 ÷ 5
            while m < N5
                for i = 0:m-1
                    j1 = id + 5k * i
                    j2 = j1 + k - 1
                    ψ1, ψ2, ψ3, ψ4 = Ψ5[m+i+1], Ψ5[2m+i+1], Ψ5[3m+i+1], Ψ5[4m+i+1]
                    for j = j1:j2
                        t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                        t1 = lazy_Mmul(a[j+k], ψ1, Q)
                        t2 = lazy_Mmul(a[j+2k], ψ2, Q)
                        t3 = lazy_Mmul(a[j+3k], ψ3, Q)
                        t4 = lazy_Mmul(a[j+4k], ψ4, Q)

                        x0 = t0 + t1
                        x0 ≥ twoQ && (x0 -= twoQ)
                        x0 = x0 + t2
                        x0 ≥ twoQ && (x0 -= twoQ)
                        x0 = x0 + t3
                        x0 ≥ twoQ && (x0 -= twoQ)
                        a[j] = x0 + t4
                        t4 = twoQ - t4

                        t0 += t4
                        t1 += t4
                        t2 += t4
                        t3 += t4

                        t0 ≥ twoQ && (t0 -= twoQ)
                        t1 ≥ twoQ && (t1 -= twoQ)
                        t2 ≥ twoQ && (t2 -= twoQ)
                        t3 ≥ twoQ && (t3 -= twoQ)

                        x1 = t0 + lazy_Mmul(t1, ω, Q)
                        x2 = t0 + lazy_Mmul(t1, ω2, Q)
                        x3 = t0 + lazy_Mmul(t1, ω3, Q)
                        x4 = t0 + lazy_Mmul(t1, ω4, Q)

                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        x1 += lazy_Mmul(t2, ω2, Q)
                        x2 += lazy_Mmul(t2, ω4, Q)
                        x3 += lazy_Mmul(t2, ω, Q)
                        x4 += lazy_Mmul(t2, ω3, Q)

                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        x1 += lazy_Mmul(t3, ω3, Q)
                        x2 += lazy_Mmul(t3, ω, Q)
                        x3 += lazy_Mmul(t3, ω4, Q)
                        x4 += lazy_Mmul(t3, ω2, Q)

                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        a[j+k] = x1
                        a[j+2k] = x2
                        a[j+3k] = x3
                        a[j+4k] = x4
                    end
                end
                m *= 5
                k ÷= 5
            end
        end
    end

    if N7 > 1
        ω, ω2, ω3, ω4, ω5, ω6 = Ψ7[end-5], Ψ7[end-4], Ψ7[end-3], Ψ7[end-2], Ψ7[end-1], Ψ7[end]
        m, k = 1, N ÷ 7
        while m < N7
            @inbounds for i = 0:m-1
                j1 = 7k * i + 1
                j2 = j1 + k - 1
                ψ1, ψ2, ψ3, ψ4, ψ5, ψ6 = Ψ7[m+i+1], Ψ7[2m+i+1], Ψ7[3m+i+1], Ψ7[4m+i+1], Ψ7[5m+i+1], Ψ7[6m+i+1]
                for j = j1:j2
                    t0 = a[j] ≥ twoQ ? a[j] - twoQ : a[j]
                    t1 = lazy_Mmul(a[j+k], ψ1, Q)
                    t2 = lazy_Mmul(a[j+2k], ψ2, Q)
                    t3 = lazy_Mmul(a[j+3k], ψ3, Q)
                    t4 = lazy_Mmul(a[j+4k], ψ4, Q)
                    t5 = lazy_Mmul(a[j+5k], ψ5, Q)
                    t6 = lazy_Mmul(a[j+6k], ψ6, Q)

                    x0 = t0 + t1
                    x0 ≥ twoQ && (x0 -= twoQ)
                    x0 = x0 + t2
                    x0 ≥ twoQ && (x0 -= twoQ)
                    x0 = x0 + t3
                    x0 ≥ twoQ && (x0 -= twoQ)
                    x0 = x0 + t4
                    x0 ≥ twoQ && (x0 -= twoQ)
                    x0 = x0 + t5
                    x0 ≥ twoQ && (x0 -= twoQ)
                    a[j] = x0 + t6
                    t6 = twoQ - t6

                    t0 += t6
                    t1 += t6
                    t2 += t6
                    t3 += t6
                    t4 += t6
                    t5 += t6

                    t0 ≥ twoQ && (t0 -= twoQ)
                    t1 ≥ twoQ && (t1 -= twoQ)
                    t2 ≥ twoQ && (t2 -= twoQ)
                    t3 ≥ twoQ && (t3 -= twoQ)
                    t4 ≥ twoQ && (t4 -= twoQ)
                    t5 ≥ twoQ && (t5 -= twoQ)

                    x1 = t0 + lazy_Mmul(t1, ω, Q)
                    x2 = t0 + lazy_Mmul(t1, ω2, Q)
                    x3 = t0 + lazy_Mmul(t1, ω3, Q)
                    x4 = t0 + lazy_Mmul(t1, ω4, Q)
                    x5 = t0 + lazy_Mmul(t1, ω5, Q)
                    x6 = t0 + lazy_Mmul(t1, ω6, Q)

                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x1 += lazy_Mmul(t2, ω2, Q)
                    x2 += lazy_Mmul(t2, ω4, Q)
                    x3 += lazy_Mmul(t2, ω6, Q)
                    x4 += lazy_Mmul(t2, ω, Q)
                    x5 += lazy_Mmul(t2, ω3, Q)
                    x6 += lazy_Mmul(t2, ω5, Q)

                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x1 += lazy_Mmul(t3, ω3, Q)
                    x2 += lazy_Mmul(t3, ω6, Q)
                    x3 += lazy_Mmul(t3, ω2, Q)
                    x4 += lazy_Mmul(t3, ω5, Q)
                    x5 += lazy_Mmul(t3, ω, Q)
                    x6 += lazy_Mmul(t3, ω4, Q)

                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x1 += lazy_Mmul(t4, ω4, Q)
                    x2 += lazy_Mmul(t4, ω, Q)
                    x3 += lazy_Mmul(t4, ω5, Q)
                    x4 += lazy_Mmul(t4, ω2, Q)
                    x5 += lazy_Mmul(t4, ω6, Q)
                    x6 += lazy_Mmul(t4, ω3, Q)

                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x1 += lazy_Mmul(t5, ω5, Q)
                    x2 += lazy_Mmul(t5, ω3, Q)
                    x3 += lazy_Mmul(t5, ω, Q)
                    x4 += lazy_Mmul(t5, ω6, Q)
                    x5 += lazy_Mmul(t5, ω4, Q)
                    x6 += lazy_Mmul(t5, ω2, Q)

                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    a[j+k] = x1
                    a[j+2k] = x2
                    a[j+3k] = x3
                    a[j+4k] = x4
                    a[j+5k] = x5
                    a[j+6k] = x6
                end
            end
            m *= 7
            k ÷= 7
        end
    end

    iMform!(a, Q)

    return nothing
end

"""
An optimised normal input Gentleman-Sande algorithm.
"""
function _intt_2a3b5c7d!(a::AbstractVector{UInt64}, Ψ2inv::Vector{UInt64}, Ψ3inv::Vector{UInt64}, Ψ5inv::Vector{UInt64}, Ψ7inv::Vector{UInt64}, Q::Modulus)::Nothing
    N = length(a)
    N2, N3, N5, N7 = factor2357(N)
    N23, N235 = N2 * N3, N2 * N3 * N5

    Mform!(a, Q)
    twoQ = Q.Q << 1
    if N2 > 1
        @inbounds for id = 1:N2:N
            m, logkp1, k = N2 >> 1, 1, 1
            for i = 0:m-1
                j = i << 1 + id
                t0 = a[j]
                u0 = a[j+k]

                a[j] = t0 + u0
                a[j] ≥ twoQ && (a[j] -= twoQ)
                a[j+k] = lazy_Mmul(t0 - u0 + twoQ, Ψ2inv[m+i+1], Q)
            end
            m >>= 1
            logkp1 += 1
            k <<= 1

            if m > 0
                for i = 0:m-1
                    j = i << 2 + id

                    t0 = a[j]
                    t1 = a[j+1]

                    u0 = a[j+k]
                    u1 = a[j+k+1]

                    a[j] = t0 + u0
                    a[j+1] = t1 + u1
                    a[j] ≥ twoQ && (a[j] -= twoQ)
                    a[j+1] ≥ twoQ && (a[j+1] -= twoQ)

                    a[j+k] = lazy_Mmul(t0 - u0 + twoQ, Ψ2inv[m+i+1], Q)
                    a[j+k+1] = lazy_Mmul(t1 - u1 + twoQ, Ψ2inv[m+i+1], Q)
                end
                m >>= 1
                logkp1 += 1
                k <<= 1
            end

            if m > 0
                for i = 0:m-1
                    j = i << 3 + id

                    t0 = a[j]
                    t1 = a[j+1]
                    t2 = a[j+2]
                    t3 = a[j+3]

                    u0 = a[j+k]
                    u1 = a[j+k+1]
                    u2 = a[j+k+2]
                    u3 = a[j+k+3]

                    a[j] = t0 + u0
                    a[j+1] = t1 + u1
                    a[j+2] = t2 + u2
                    a[j+3] = t3 + u3
                    a[j] ≥ twoQ && (a[j] -= twoQ)
                    a[j+1] ≥ twoQ && (a[j+1] -= twoQ)
                    a[j+2] ≥ twoQ && (a[j+2] -= twoQ)
                    a[j+3] ≥ twoQ && (a[j+3] -= twoQ)

                    a[j+k] = lazy_Mmul(t0 - u0 + twoQ, Ψ2inv[m+i+1], Q)
                    a[j+k+1] = lazy_Mmul(t1 - u1 + twoQ, Ψ2inv[m+i+1], Q)
                    a[j+k+2] = lazy_Mmul(t2 - u2 + twoQ, Ψ2inv[m+i+1], Q)
                    a[j+k+3] = lazy_Mmul(t3 - u3 + twoQ, Ψ2inv[m+i+1], Q)
                end
                m >>= 1
                logkp1 += 1
                k <<= 1
            end

            while m > 1
                for i = 0:m-1
                    j1 = i << logkp1 + id
                    j2 = j1 + k - 1
                    for j = j1:8:j2
                        t0 = a[j]
                        t1 = a[j+1]
                        t2 = a[j+2]
                        t3 = a[j+3]
                        t4 = a[j+4]
                        t5 = a[j+5]
                        t6 = a[j+6]
                        t7 = a[j+7]

                        u0 = a[j+k]
                        u1 = a[j+k+1]
                        u2 = a[j+k+2]
                        u3 = a[j+k+3]
                        u4 = a[j+k+4]
                        u5 = a[j+k+5]
                        u6 = a[j+k+6]
                        u7 = a[j+k+7]

                        a[j] = t0 + u0
                        a[j+1] = t1 + u1
                        a[j+2] = t2 + u2
                        a[j+3] = t3 + u3
                        a[j+4] = t4 + u4
                        a[j+5] = t5 + u5
                        a[j+6] = t6 + u6
                        a[j+7] = t7 + u7
                        a[j] ≥ twoQ ? (a[j] -= twoQ) : a[j]
                        a[j+1] ≥ twoQ ? (a[j+1] -= twoQ) : a[j+1]
                        a[j+2] ≥ twoQ ? (a[j+2] -= twoQ) : a[j+2]
                        a[j+3] ≥ twoQ ? (a[j+3] -= twoQ) : a[j+3]
                        a[j+4] ≥ twoQ ? (a[j+4] -= twoQ) : a[j+4]
                        a[j+5] ≥ twoQ ? (a[j+5] -= twoQ) : a[j+5]
                        a[j+6] ≥ twoQ ? (a[j+6] -= twoQ) : a[j+6]
                        a[j+7] ≥ twoQ ? (a[j+7] -= twoQ) : a[j+7]

                        a[j+k] = lazy_Mmul(t0 - u0 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+1] = lazy_Mmul(t1 - u1 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+2] = lazy_Mmul(t2 - u2 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+3] = lazy_Mmul(t3 - u3 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+4] = lazy_Mmul(t4 - u4 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+5] = lazy_Mmul(t5 - u5 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+6] = lazy_Mmul(t6 - u6 + twoQ, Ψ2inv[m+i+1], Q)
                        a[j+k+7] = lazy_Mmul(t7 - u7 + twoQ, Ψ2inv[m+i+1], Q)
                    end
                end
                m >>= 1
                logkp1 += 1
                k <<= 1
            end

            if m > 0
                for j = id:id+k-1
                    t0 = a[j]
                    u0 = a[j+k]
                    a[j] = t0 + u0
                    a[j] ≥ twoQ && (a[j] -= twoQ)
                    a[j+k] = lazy_Mmul(t0 - u0 + twoQ, Ψ2inv[2], Q)
                end
            end
        end
    end

    if N3 > 1
        ω, ω2 = Ψ3inv[end-1], Ψ3inv[end]
        for id = 1:N23:N
            m, k = N3 ÷ 3, N2
            while m > 0
                for i = 0:m-1
                    j1 = id + 3k * i
                    j2 = j1 + k - 1
                    ψ1, ψ2 = Ψ3inv[m+i+1], Ψ3inv[2m+i+1]
                    for j = j1:j2
                        t0 = a[j]
                        t1 = a[j+k]
                        t2 = a[j+2k]

                        x0 = t0 + t1
                        x1 = t0 + lazy_Mmul(t1, ω2, Q)
                        x2 = t0 + lazy_Mmul(t1, ω, Q)

                        x0 ≥ twoQ && (x0 -= twoQ)
                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)

                        x0 += t2
                        x1 += lazy_Mmul(t2, ω, Q)
                        x2 += lazy_Mmul(t2, ω2, Q)

                        x0 ≥ twoQ && (x0 -= twoQ)
                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)

                        a[j] = x0
                        a[j+k] = lazy_Mmul(x1, ψ1, Q)
                        a[j+2k] = lazy_Mmul(x2, ψ2, Q)
                    end
                end
                m ÷= 3
                k *= 3
            end
        end
    end

    if N5 > 1
        ω, ω2, ω3, ω4 = Ψ5inv[end-3], Ψ5inv[end-2], Ψ5inv[end-1], Ψ5inv[end]
        for id = 1:N235:N
            m, k = N5 ÷ 5, N23
            while m > 0
                for i = 0:m-1
                    j1 = id + 5k * i
                    j2 = j1 + k - 1
                    ψ1, ψ2, ψ3, ψ4 = Ψ5inv[m+i+1], Ψ5inv[2m+i+1], Ψ5inv[3m+i+1], Ψ5inv[4m+i+1]
                    for j = j1:j2
                        t0 = a[j]
                        t1 = a[j+k]
                        t2 = a[j+2k]
                        t3 = a[j+3k]
                        t4 = a[j+4k]

                        x0 = t0 + t1
                        x1 = t0 + lazy_Mmul(t1, ω4, Q)
                        x2 = t0 + lazy_Mmul(t1, ω3, Q)
                        x3 = t0 + lazy_Mmul(t1, ω2, Q)
                        x4 = t0 + lazy_Mmul(t1, ω, Q)

                        x0 ≥ twoQ && (x0 -= twoQ)
                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        x0 += t2
                        x1 += lazy_Mmul(t2, ω3, Q)
                        x2 += lazy_Mmul(t2, ω, Q)
                        x3 += lazy_Mmul(t2, ω4, Q)
                        x4 += lazy_Mmul(t2, ω2, Q)

                        x0 ≥ twoQ && (x0 -= twoQ)
                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        x0 += t3
                        x1 += lazy_Mmul(t3, ω2, Q)
                        x2 += lazy_Mmul(t3, ω4, Q)
                        x3 += lazy_Mmul(t3, ω, Q)
                        x4 += lazy_Mmul(t3, ω3, Q)

                        x0 ≥ twoQ && (x0 -= twoQ)
                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        x0 += t4
                        x1 += lazy_Mmul(t4, ω, Q)
                        x2 += lazy_Mmul(t4, ω2, Q)
                        x3 += lazy_Mmul(t4, ω3, Q)
                        x4 += lazy_Mmul(t4, ω4, Q)

                        x0 ≥ twoQ && (x0 -= twoQ)
                        x1 ≥ twoQ && (x1 -= twoQ)
                        x2 ≥ twoQ && (x2 -= twoQ)
                        x3 ≥ twoQ && (x3 -= twoQ)
                        x4 ≥ twoQ && (x4 -= twoQ)

                        a[j] = x0
                        a[j+k] = lazy_Mmul(x1, ψ1, Q)
                        a[j+2k] = lazy_Mmul(x2, ψ2, Q)
                        a[j+3k] = lazy_Mmul(x3, ψ3, Q)
                        a[j+4k] = lazy_Mmul(x4, ψ4, Q)
                    end
                end
                m ÷= 5
                k *= 5
            end
        end
    end

    if N7 > 1
        ω, ω2, ω3, ω4, ω5, ω6 = Ψ7inv[end-5], Ψ7inv[end-4], Ψ7inv[end-3], Ψ7inv[end-2], Ψ7inv[end-1], Ψ7inv[end]
        m, k = N7 ÷ 7, N235
        while m > 0
            for i = 0:m-1
                j1 = 7k * i + 1
                j2 = j1 + k - 1
                ψ1, ψ2, ψ3, ψ4, ψ5, ψ6 = Ψ7inv[m+i+1], Ψ7inv[2m+i+1], Ψ7inv[3m+i+1], Ψ7inv[4m+i+1], Ψ7inv[5m+i+1], Ψ7inv[6m+i+1]
                for j = j1:j2
                    t0 = a[j]
                    t1 = a[j+k]
                    t2 = a[j+2k]
                    t3 = a[j+3k]
                    t4 = a[j+4k]
                    t5 = a[j+5k]
                    t6 = a[j+6k]

                    x0 = t0 + t1
                    x1 = t0 + lazy_Mmul(t1, ω6, Q)
                    x2 = t0 + lazy_Mmul(t1, ω5, Q)
                    x3 = t0 + lazy_Mmul(t1, ω4, Q)
                    x4 = t0 + lazy_Mmul(t1, ω3, Q)
                    x5 = t0 + lazy_Mmul(t1, ω2, Q)
                    x6 = t0 + lazy_Mmul(t1, ω, Q)

                    x0 ≥ twoQ && (x0 -= twoQ)
                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x0 += t2
                    x1 += lazy_Mmul(t2, ω5, Q)
                    x2 += lazy_Mmul(t2, ω3, Q)
                    x3 += lazy_Mmul(t2, ω, Q)
                    x4 += lazy_Mmul(t2, ω6, Q)
                    x5 += lazy_Mmul(t2, ω4, Q)
                    x6 += lazy_Mmul(t2, ω2, Q)

                    x0 ≥ twoQ && (x0 -= twoQ)
                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x0 += t3
                    x1 += lazy_Mmul(t3, ω4, Q)
                    x2 += lazy_Mmul(t3, ω, Q)
                    x3 += lazy_Mmul(t3, ω5, Q)
                    x4 += lazy_Mmul(t3, ω2, Q)
                    x5 += lazy_Mmul(t3, ω6, Q)
                    x6 += lazy_Mmul(t3, ω3, Q)

                    x0 ≥ twoQ && (x0 -= twoQ)
                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x0 += t4
                    x1 += lazy_Mmul(t4, ω3, Q)
                    x2 += lazy_Mmul(t4, ω6, Q)
                    x3 += lazy_Mmul(t4, ω2, Q)
                    x4 += lazy_Mmul(t4, ω5, Q)
                    x5 += lazy_Mmul(t4, ω, Q)
                    x6 += lazy_Mmul(t4, ω4, Q)

                    x0 ≥ twoQ && (x0 -= twoQ)
                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x0 += t5
                    x1 += lazy_Mmul(t5, ω2, Q)
                    x2 += lazy_Mmul(t5, ω4, Q)
                    x3 += lazy_Mmul(t5, ω6, Q)
                    x4 += lazy_Mmul(t5, ω, Q)
                    x5 += lazy_Mmul(t5, ω3, Q)
                    x6 += lazy_Mmul(t5, ω5, Q)

                    x0 ≥ twoQ && (x0 -= twoQ)
                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    x0 += t6
                    x1 += lazy_Mmul(t6, ω, Q)
                    x2 += lazy_Mmul(t6, ω2, Q)
                    x3 += lazy_Mmul(t6, ω3, Q)
                    x4 += lazy_Mmul(t6, ω4, Q)
                    x5 += lazy_Mmul(t6, ω5, Q)
                    x6 += lazy_Mmul(t6, ω6, Q)

                    x0 ≥ twoQ && (x0 -= twoQ)
                    x1 ≥ twoQ && (x1 -= twoQ)
                    x2 ≥ twoQ && (x2 -= twoQ)
                    x3 ≥ twoQ && (x3 -= twoQ)
                    x4 ≥ twoQ && (x4 -= twoQ)
                    x5 ≥ twoQ && (x5 -= twoQ)
                    x6 ≥ twoQ && (x6 -= twoQ)

                    a[j] = x0
                    a[j+k] = lazy_Mmul(x1, ψ1, Q)
                    a[j+2k] = lazy_Mmul(x2, ψ2, Q)
                    a[j+3k] = lazy_Mmul(x3, ψ3, Q)
                    a[j+4k] = lazy_Mmul(x4, ψ4, Q)
                    a[j+5k] = lazy_Mmul(x5, ψ5, Q)
                    a[j+6k] = lazy_Mmul(x6, ψ6, Q)
                end
            end
            m ÷= 7
            k *= 7
        end
    end

    iMform!(a, Q)

    return nothing
end

struct CyclicNTTransformer2a3b5c7d <: CyclicNTTransformer
    Q::Modulus
    N::Int64
    N⁻¹::UInt64
    Ψ2::Vector{UInt64}
    Ψ3::Vector{UInt64}
    Ψ5::Vector{UInt64}
    Ψ7::Vector{UInt64}
    Ψ2inv::Vector{UInt64}
    Ψ3inv::Vector{UInt64}
    Ψ5inv::Vector{UInt64}
    Ψ7inv::Vector{UInt64}
    iseasy::Bool
    idx::Vector{Int64}
    buff::Vector{UInt64}

    function CyclicNTTransformer2a3b5c7d(N::Int64, Q::Modulus)::CyclicNTTransformer2a3b5c7d
        if Q.Q % N ≠ 1
            throw(DomainError("Modulus does not have enough factors to perform NTT."))
        end

        N2, N3, N5, N7 = factor2357(N)

        ξ = primitive_root_finder(Q)

        Ψ2, Ψ2inv = gen_roots(N2, 2, ξ, Q)
        Ψ3, Ψ3inv = gen_roots(N3, 3, ξ, Q)
        Ψ5, Ψ5inv = gen_roots(N5, 5, ξ, Q)
        Ψ7, Ψ7inv = gen_roots(N7, 7, ξ, Q)

        iseasy = (N2 == N) || (N3 == N) || (N5 == N) || (N7 == N)

        if !iseasy
            idx = Vector{Int64}(undef, N)
            @inbounds for i1 = 0:N2-1, i2 = 0:N3-1, i3 = 0:N5-1, i4 = 0:N7-1
                idx[i1+i2*N2+i3*N2*N3+i4*N2*N3*N5+1] = (i1 * N3 * N5 * N7 + i2 * N2 * N5 * N7 + i3 * N2 * N3 * N7 + i4 * N2 * N3 * N5) % N + 1
            end

            buff = Vector{UInt64}(undef, N)
        else
            idx = Int64[]
            buff = UInt64[]
        end

        new(Q, N, invmod(N, Q), Ψ2, Ψ3, Ψ5, Ψ7, Ψ2inv, Ψ3inv, Ψ5inv, Ψ7inv, iseasy, idx, buff)
    end
end

@views ntt!(a::AbstractVector{UInt64}, ntter::CyclicNTTransformer2a3b5c7d)::Nothing = begin
    if !ntter.iseasy
        @inbounds @simd for i = eachindex(a)
            ntter.buff[i] = a[ntter.idx[i]]
        end
        @inbounds @simd for i = eachindex(a)
            a[i] = ntter.buff[i]
        end
    end

    _ntt_2a3b5c7d!(a, ntter.Ψ2, ntter.Ψ3, ntter.Ψ5, ntter.Ψ7, ntter.Q)

    return nothing
end

@views ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclicNTTransformer2a3b5c7d)::Nothing = begin
    @. res = a
    ntt!(res, ntter)

    return nothing
end

@views intt!(a::AbstractVector{UInt64}, ntter::CyclicNTTransformer2a3b5c7d)::Nothing = begin
    _intt_2a3b5c7d!(a, ntter.Ψ2inv, ntter.Ψ3inv, ntter.Ψ5inv, ntter.Ψ7inv, ntter.Q)

    if !ntter.iseasy
        @inbounds @simd for i = eachindex(a)
            ntter.buff[ntter.idx[i]] = a[i]
        end
        @inbounds @simd for i = eachindex(a)
            a[i] = ntter.buff[i]
        end
    end

    Bmul_to!(a, ntter.N⁻¹, a, ntter.Q)

    return nothing
end

@views intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclicNTTransformer2a3b5c7d)::Nothing = begin
    @. res = a
    intt!(res, ntter)

    return nothing
end

struct CyclicNTTransformerBluestein <: CyclicNTTransformer      
    Q::Modulus
    m::Int64
    N::Int64
    Ψ2::Vector{UInt64}
    Ψ3::Vector{UInt64}
    Ψ5::Vector{UInt64}
    Ψ7::Vector{UInt64}
    Ψ2inv::Vector{UInt64}
    Ψ3inv::Vector{UInt64}
    Ψ5inv::Vector{UInt64}
    Ψ7inv::Vector{UInt64}
    iseasy::Bool
    idx::Vector{Int64}
    buff::Vector{UInt64}
    Ζ::Vector{UInt64}
    Ζinv::Vector{UInt64}
    chirp::Vector{UInt64}
    chirpinv::Vector{UInt64}

    function CyclicNTTransformerBluestein(m::Int64, Q::Modulus)::CyclicNTTransformerBluestein
        N = next2a3b(2m - 1)
        N2, N3, N5, N7 = factor2357(N)

        if Q.Q % 2m ≠ 1
            throw(DomainError("Modulus does not have $(2m)-th root of unity."))
        end
        if Q.Q % N ≠ 1
            throw(DomainError("Modulus does not have enough two factor for Bluestein NTT."))
        end

        ξ = primitive_root_finder(Q)

        Ψ2, Ψ2inv = gen_roots(N2, 2, ξ, Q)
        Ψ3, Ψ3inv = gen_roots(N3, 3, ξ, Q)
        Ψ5, Ψ5inv = gen_roots(N5, 5, ξ, Q)
        Ψ7, Ψ7inv = gen_roots(N7, 7, ξ, Q)

        idx = Int64[i^2 % 2m for i = 0:m-1]
        ζ = ith_root_from_primitive_root(2m, ξ, Q)
        ζinv = powermod(ζ, 2m - 1, Q)
        Ζ, Ζinv = powermod.(ζ, idx, Ref(Q)), powermod.(ζinv, idx, Ref(Q))

        iseasy = (N2 == N) || (N3 == N) || (N5 == N) || (N7 == N)

        if !iseasy
            idx = Vector{Int64}(undef, N)
            for i1 = 0:N2-1, i2 = 0:N3-1, i3 = 0:N5-1, i4 = 0:N7-1
                idx[i1+i2*N2+i3*N2*N3+i4*N2*N3*N5+1] = (i1 * N3 * N5 * N7 + i2 * N2 * N5 * N7 + i3 * N2 * N3 * N7 + i4 * N2 * N3 * N5) % N + 1
            end
        else
            idx = Int64[]
        end

        buff = Vector{UInt64}(undef, N)

        # We scale chirp and chirpinv instead of dividing by N for convolution, for simplicity.
        chirp = vcat(Ζinv, zeros(UInt64, N - 2m + 1), Ζinv[end:-1:2])
        Bmul_to!(chirp, invmod(N, Q), chirp, Q)
        chirpinv = vcat(Ζ, zeros(UInt64, N - 2m + 1), Ζ[end:-1:2])
        Bmul_to!(chirpinv, invmod(N * m, Q), chirpinv, Q)

        if !iseasy
            @inbounds @simd for i = eachindex(chirp)
                buff[i] = chirp[idx[i]]
            end
            @. chirp = buff
        end
        _ntt_2a3b5c7d!(chirp, Ψ2, Ψ3, Ψ5, Ψ7, Q)

        if !iseasy
            @inbounds @simd for i = eachindex(chirpinv)
                buff[i] = chirpinv[idx[i]]
            end
            @. chirpinv = buff
        end
        _ntt_2a3b5c7d!(chirpinv, Ψ2, Ψ3, Ψ5, Ψ7, Q)

        new(Q, m, N, Ψ2, Ψ3, Ψ5, Ψ7, Ψ2inv, Ψ3inv, Ψ5inv, Ψ7inv, iseasy, idx, buff, Ζ, Ζinv, chirp, chirpinv)
    end
end

# Bluestein NTT
@views function ntt!(a::AbstractVector{UInt64}, ntter::CyclicNTTransformerBluestein)::Nothing
    @. ntter.buff = 0
    if ntter.iseasy
        @inbounds for i = 1:ntter.m
            ntter.buff[i] = lazy_Bmul(a[i], ntter.Ζ[i], ntter.Q)
        end
    else
        @inbounds for i = 1:ntter.N
            ntter.idx[i] ≤ ntter.m && (ntter.buff[i] = lazy_Bmul(a[ntter.idx[i]], ntter.Ζ[ntter.idx[i]], ntter.Q))
        end
    end

    _ntt_2a3b5c7d!(ntter.buff, ntter.Ψ2, ntter.Ψ3, ntter.Ψ5, ntter.Ψ7, ntter.Q)
    lazy_Bmul_to!(ntter.buff, ntter.buff, ntter.chirp, ntter.Q)
    _intt_2a3b5c7d!(ntter.buff, ntter.Ψ2inv, ntter.Ψ3inv, ntter.Ψ5inv, ntter.Ψ7inv, ntter.Q)

    if ntter.iseasy
        @inbounds for i = 1:ntter.m
            a[i] = Bmul(ntter.buff[i], ntter.Ζ[i], ntter.Q)
        end
    else
        @inbounds for i = 1:ntter.N
            ntter.idx[i] ≤ ntter.m && (a[ntter.idx[i]] = Bmul(ntter.buff[i], ntter.Ζ[ntter.idx[i]], ntter.Q))
        end
    end

    return nothing
end

@views ntt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclicNTTransformerBluestein)::Nothing = begin
    @. res = a
    ntt!(res, ntter)

    return nothing
end

# Bluestein iNTT
@views function intt!(a::AbstractVector{UInt64}, ntter::CyclicNTTransformerBluestein)::Nothing
    @. ntter.buff = 0
    if ntter.iseasy
        @inbounds for i = 1:ntter.m
            ntter.buff[i] = lazy_Bmul(a[i], ntter.Ζinv[i], ntter.Q)
        end
    else
        @inbounds for i = 1:ntter.N
            ntter.idx[i] ≤ ntter.m && (ntter.buff[i] = lazy_Bmul(a[ntter.idx[i]], ntter.Ζinv[ntter.idx[i]], ntter.Q))
        end
    end

    _ntt_2a3b5c7d!(ntter.buff, ntter.Ψ2, ntter.Ψ3, ntter.Ψ5, ntter.Ψ7, ntter.Q)
    lazy_Bmul_to!(ntter.buff, ntter.buff, ntter.chirpinv, ntter.Q)
    _intt_2a3b5c7d!(ntter.buff, ntter.Ψ2inv, ntter.Ψ3inv, ntter.Ψ5inv, ntter.Ψ7inv, ntter.Q)

    if ntter.iseasy
        @inbounds for i = 1:ntter.m
            a[i] = Bmul(ntter.buff[i], ntter.Ζinv[i], ntter.Q)
        end
    else
        @inbounds for i = 1:ntter.N
            ntter.idx[i] ≤ ntter.m && (a[ntter.idx[i]] = Bmul(ntter.buff[i], ntter.Ζinv[ntter.idx[i]], ntter.Q))
        end
    end

    return nothing
end

@views intt_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ntter::CyclicNTTransformerBluestein)::Nothing = begin
    @. res = a
    intt!(res, ntter)

    return nothing
end

(::Type{CyclicNTTransformer})(m::Int64, Q::Modulus)::CyclicNTTransformer =
    is2a3b5c7d(m) ? CyclicNTTransformer2a3b5c7d(m, Q) : CyclicNTTransformerBluestein(m, Q)