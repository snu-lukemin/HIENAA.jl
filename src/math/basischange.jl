"""
BasisExtender provides basis extension from Q to P.
"""
struct BasisExtender
    Q::Vector{Modulus}
    P::Vector{Modulus}
    Qtilde::Vector{UInt64}
    Qstar::Array{UInt64,2}
    QmodP::Vector{UInt64}
    oneoverQ::Vector{UInt128}
    isPinQ::Vector{Int64}
    v128::Vector{UInt128}
    v64::Vector{UInt64}
    buff::Array{UInt64,2}

    function BasisExtender(Q::Moduli, P::Moduli)::BasisExtender
        Qlen, Plen = length(Q), length(P)

        Qtilde = Vector{UInt64}(undef, Qlen)
        Qstar = Array{UInt64,2}(undef, Qlen, Plen)
        QmodP = Vector{UInt64}(undef, Plen)
        oneoverQ = Vector{UInt128}(undef, Qlen)
        isPinQ = zeros(Int64, Plen)
        v128 = Vector{UInt128}(undef, 8)
        v64 = Vector{UInt64}(undef, 8)
        buff = Array{UInt64}(undef, 8, Qlen)

        setprecision(192)
        @inbounds for i = 1:Qlen
            Qtilde[i] = 1
            for j = 1:Qlen
                i == j && continue
                Qtilde[i] = Bmul(Qtilde[i], invmod(Q[j].Q, Q[i]), Q[i])
            end

            oneoverQ[i] = round(UInt128, (Int128(1) << fixed_prec) / big(Q[i].Q))
        end

        @inbounds for j = 1:Plen
            QmodP[j] = 1

            for i = 1:Qlen
                Qstar[i, j] = 1
                for k = 1:Qlen
                    i == k && continue
                    Qstar[i, j] = Bmul(Qstar[i, j], Q[k].Q, P[j])
                end

                QmodP[j] = Bmul(QmodP[j], Q[i].Q, P[j])

                if P[j].Q == Q[i].Q
                    isPinQ[j] = i
                    break
                end
            end
        end

        new(collect(Q), collect(P), Qtilde, Qstar, QmodP, oneoverQ, isPinQ, v128, v64, buff)
    end
end

"""
basis_extend_to! performs the HPS18 Basis Extension algorithm.
"""
function basis_extend_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, be::BasisExtender)::Nothing
    Q, P = be.Q, be.P
    Qlen, Plen = min(length(a), length(Q)), min(length(res), length(P))

    Qtilde, Qstar, QmodP, oneoverQ, isPinQ, buff = be.Qtilde, be.Qstar, be.QmodP, be.oneoverQ, be.isPinQ, be.buff

    @inbounds for i = 1:Qlen
        buff[i] = Bmul(a[i], Qtilde[i], Q[i])
    end

    v128 = zero(UInt128)
    @inbounds for i = 1:Qlen
        v128 += mult_and_round(buff[i], oneoverQ[i])
    end
    v64 = round_to_uint64(v128)

    @inbounds for j = 1:Plen
        # Decompose more efficiently.
        if isPinQ[j] > 0
            res[j] = a[isPinQ[j]]
            continue
        end

        res[j] = lazy_Bmul(neg(QmodP[j], P[j]), v64, P[j])
        for i = 1:Qlen-1
            res[j] = lazy_Bred(widemul(buff[i], Qstar[i, j]) + res[j], P[j])
        end
        res[j] = Bred(widemul(buff[Qlen], Qstar[Qlen, j]) + res[j], P[j])
    end

    return nothing
end

Base.@propagate_inbounds function basis_extend_to!(res::AbstractVector{Vector{UInt64}}, a::AbstractVector{Vector{UInt64}}, be::BasisExtender)::Nothing
    N, reslen, alen = length(res[1]), length(res), length(a)

    Q, P = be.Q, be.P
    Qlen, Plen = min(alen, length(Q)), min(reslen, length(P))

    if Qlen == 1
        for i = 1:Plen
            if be.isPinQ[i] == 1
                copy!(res[i], a[1])
            else
                halfQ = Q[1].Q >> 1
                for k = 0:8:N-8
                    res[i][k+1] = a[1][k+1] ≤ halfQ ? Bred(a[1][k+1], P[i]) : Bred(signed(a[1][k+1] - Q[1].Q), P[i])
                    res[i][k+2] = a[1][k+2] ≤ halfQ ? Bred(a[1][k+2], P[i]) : Bred(signed(a[1][k+2] - Q[1].Q), P[i])
                    res[i][k+3] = a[1][k+3] ≤ halfQ ? Bred(a[1][k+3], P[i]) : Bred(signed(a[1][k+3] - Q[1].Q), P[i])
                    res[i][k+4] = a[1][k+4] ≤ halfQ ? Bred(a[1][k+4], P[i]) : Bred(signed(a[1][k+4] - Q[1].Q), P[i])
                    res[i][k+5] = a[1][k+5] ≤ halfQ ? Bred(a[1][k+5], P[i]) : Bred(signed(a[1][k+5] - Q[1].Q), P[i])
                    res[i][k+6] = a[1][k+6] ≤ halfQ ? Bred(a[1][k+6], P[i]) : Bred(signed(a[1][k+6] - Q[1].Q), P[i])
                    res[i][k+7] = a[1][k+7] ≤ halfQ ? Bred(a[1][k+7], P[i]) : Bred(signed(a[1][k+7] - Q[1].Q), P[i])
                    res[i][k+8] = a[1][k+8] ≤ halfQ ? Bred(a[1][k+8], P[i]) : Bred(signed(a[1][k+8] - Q[1].Q), P[i])
                end

                for k = 8(N>>3)+1:N
                    res[i][k] = a[1][k] ≤ halfQ ? Bred(a[1][k], P[i]) : Bred(signed(a[1][k] - Q[1].Q), P[i])
                end
            end
        end
    else
        Qtilde, Qstar, QmodP, oneoverQ, isPinQ, v128, v64, buff = be.Qtilde, be.Qstar, be.QmodP, be.oneoverQ, be.isPinQ, be.v128, be.v64, be.buff

        for k = 0:8:N-8
            for i = 1:Qlen
                buff[1, i] = Bmul(a[i][k+1], Qtilde[i], Q[i])
                buff[2, i] = Bmul(a[i][k+2], Qtilde[i], Q[i])
                buff[3, i] = Bmul(a[i][k+3], Qtilde[i], Q[i])
                buff[4, i] = Bmul(a[i][k+4], Qtilde[i], Q[i])
                buff[5, i] = Bmul(a[i][k+5], Qtilde[i], Q[i])
                buff[6, i] = Bmul(a[i][k+6], Qtilde[i], Q[i])
                buff[7, i] = Bmul(a[i][k+7], Qtilde[i], Q[i])
                buff[8, i] = Bmul(a[i][k+8], Qtilde[i], Q[i])
            end

            v128[1] = 0
            v128[2] = 0
            v128[3] = 0
            v128[4] = 0
            v128[5] = 0
            v128[6] = 0
            v128[7] = 0
            v128[8] = 0
            for i = 1:Qlen
                v128[1] += mult_and_round(buff[1, i], oneoverQ[i])
                v128[2] += mult_and_round(buff[2, i], oneoverQ[i])
                v128[3] += mult_and_round(buff[3, i], oneoverQ[i])
                v128[4] += mult_and_round(buff[4, i], oneoverQ[i])
                v128[5] += mult_and_round(buff[5, i], oneoverQ[i])
                v128[6] += mult_and_round(buff[6, i], oneoverQ[i])
                v128[7] += mult_and_round(buff[7, i], oneoverQ[i])
                v128[8] += mult_and_round(buff[8, i], oneoverQ[i])
            end
            v64[1] = round_to_uint64(v128[1])
            v64[2] = round_to_uint64(v128[2])
            v64[3] = round_to_uint64(v128[3])
            v64[4] = round_to_uint64(v128[4])
            v64[5] = round_to_uint64(v128[5])
            v64[6] = round_to_uint64(v128[6])
            v64[7] = round_to_uint64(v128[7])
            v64[8] = round_to_uint64(v128[8])

            for j = 1:Plen
                # Decompose more efficiently.
                if Qlen ≥ isPinQ[j] > 0
                    res[j][k+1] = a[isPinQ[j]][k+1]
                    res[j][k+2] = a[isPinQ[j]][k+2]
                    res[j][k+3] = a[isPinQ[j]][k+3]
                    res[j][k+4] = a[isPinQ[j]][k+4]
                    res[j][k+5] = a[isPinQ[j]][k+5]
                    res[j][k+6] = a[isPinQ[j]][k+6]
                    res[j][k+7] = a[isPinQ[j]][k+7]
                    res[j][k+8] = a[isPinQ[j]][k+8]
                    continue
                end

                negQmodP = P[j].Q - QmodP[j]
                res[j][k+1] = lazy_Bmul(negQmodP, v64[1], P[j])
                res[j][k+2] = lazy_Bmul(negQmodP, v64[2], P[j])
                res[j][k+3] = lazy_Bmul(negQmodP, v64[3], P[j])
                res[j][k+4] = lazy_Bmul(negQmodP, v64[4], P[j])
                res[j][k+5] = lazy_Bmul(negQmodP, v64[5], P[j])
                res[j][k+6] = lazy_Bmul(negQmodP, v64[6], P[j])
                res[j][k+7] = lazy_Bmul(negQmodP, v64[7], P[j])
                res[j][k+8] = lazy_Bmul(negQmodP, v64[8], P[j])
                for i = 1:Qlen-1
                    res[j][k+1] = lazy_Bred(widemul(buff[1, i], Qstar[i, j]) + res[j][k+1], P[j])
                    res[j][k+2] = lazy_Bred(widemul(buff[2, i], Qstar[i, j]) + res[j][k+2], P[j])
                    res[j][k+3] = lazy_Bred(widemul(buff[3, i], Qstar[i, j]) + res[j][k+3], P[j])
                    res[j][k+4] = lazy_Bred(widemul(buff[4, i], Qstar[i, j]) + res[j][k+4], P[j])
                    res[j][k+5] = lazy_Bred(widemul(buff[5, i], Qstar[i, j]) + res[j][k+5], P[j])
                    res[j][k+6] = lazy_Bred(widemul(buff[6, i], Qstar[i, j]) + res[j][k+6], P[j])
                    res[j][k+7] = lazy_Bred(widemul(buff[7, i], Qstar[i, j]) + res[j][k+7], P[j])
                    res[j][k+8] = lazy_Bred(widemul(buff[8, i], Qstar[i, j]) + res[j][k+8], P[j])
                end
                res[j][k+1] = Bred(widemul(buff[1, Qlen], Qstar[Qlen, j]) + res[j][k+1], P[j])
                res[j][k+2] = Bred(widemul(buff[2, Qlen], Qstar[Qlen, j]) + res[j][k+2], P[j])
                res[j][k+3] = Bred(widemul(buff[3, Qlen], Qstar[Qlen, j]) + res[j][k+3], P[j])
                res[j][k+4] = Bred(widemul(buff[4, Qlen], Qstar[Qlen, j]) + res[j][k+4], P[j])
                res[j][k+5] = Bred(widemul(buff[5, Qlen], Qstar[Qlen, j]) + res[j][k+5], P[j])
                res[j][k+6] = Bred(widemul(buff[6, Qlen], Qstar[Qlen, j]) + res[j][k+6], P[j])
                res[j][k+7] = Bred(widemul(buff[7, Qlen], Qstar[Qlen, j]) + res[j][k+7], P[j])
                res[j][k+8] = Bred(widemul(buff[8, Qlen], Qstar[Qlen, j]) + res[j][k+8], P[j])
            end
        end

        for k = 8(N>>3)+1:N
            for i = 1:Qlen
                buff[i] = Bmul(a[i][k], Qtilde[i], Q[i])
            end

            v128 = zero(UInt128)
            for i = 1:Qlen
                v128 += mult_and_round(buff[i], oneoverQ[i])
            end
            v64 = round_to_uint64(v128)

            for j = 1:Plen
                # Decompose more efficiently.
                if isPinQ[j] > 0
                    res[j][k] = a[isPinQ[j]][k]
                    continue
                end

                res[j][k] = lazy_Bmul(P[j].Q - QmodP[j], v64, P[j])
                for i = 1:Qlen-1
                    res[j][k] = lazy_Bred(widemul(buff[i], Qstar[i, j]) + res[j][k], P[j])
                end
                res[j][k] = Bred(widemul(buff[Qlen], Qstar[Qlen, j]) + res[j][k], P[j])
            end
        end
    end

    return nothing
end

function basis_extend_to!(res::AbstractVector{Vector{UInt64}}, idxres::AbstractRange{Int64}, a::AbstractVector{Vector{UInt64}}, idxa::AbstractRange{Int64}, be::BasisExtender)::Nothing
    if length(idxres) ≠ length(idxa)
        throw(DimensionMismatch("The length of the input vectors does not match."))
    end

    N, reslen, alen = length(idxres), length(res), length(a)

    Q, P = be.Q, be.P
    Qlen, Plen = min(alen, length(Q)), min(reslen, length(P))

    if Qlen == 1
        @inbounds for i = 1:Plen
            if be.isPinQ[i] == 1
                @views @. res[i][idxres] = a[1][idxa]
            else
                halfQ = Q[1].Q >> 1
                for k = 0:8:N-8
                    res[i][idxres[k+1]] = a[1][idxa[k+1]] ≤ halfQ ? Bred(a[1][idxa[k+1]], P[i]) : Bred(signed(a[1][idxa[k+1]] - Q[1].Q), P[i])
                    res[i][idxres[k+2]] = a[1][idxa[k+2]] ≤ halfQ ? Bred(a[1][idxa[k+2]], P[i]) : Bred(signed(a[1][idxa[k+2]] - Q[1].Q), P[i])
                    res[i][idxres[k+3]] = a[1][idxa[k+3]] ≤ halfQ ? Bred(a[1][idxa[k+3]], P[i]) : Bred(signed(a[1][idxa[k+3]] - Q[1].Q), P[i])
                    res[i][idxres[k+4]] = a[1][idxa[k+4]] ≤ halfQ ? Bred(a[1][idxa[k+4]], P[i]) : Bred(signed(a[1][idxa[k+4]] - Q[1].Q), P[i])
                    res[i][idxres[k+5]] = a[1][idxa[k+5]] ≤ halfQ ? Bred(a[1][idxa[k+5]], P[i]) : Bred(signed(a[1][idxa[k+5]] - Q[1].Q), P[i])
                    res[i][idxres[k+6]] = a[1][idxa[k+6]] ≤ halfQ ? Bred(a[1][idxa[k+6]], P[i]) : Bred(signed(a[1][idxa[k+6]] - Q[1].Q), P[i])
                    res[i][idxres[k+7]] = a[1][idxa[k+7]] ≤ halfQ ? Bred(a[1][idxa[k+7]], P[i]) : Bred(signed(a[1][idxa[k+7]] - Q[1].Q), P[i])
                    res[i][idxres[k+8]] = a[1][idxa[k+8]] ≤ halfQ ? Bred(a[1][idxa[k+8]], P[i]) : Bred(signed(a[1][idxa[k+8]] - Q[1].Q), P[i])
                end

                for k = 8(N>>3)+1:N
                    res[i][idxres[k]] = a[1][idxa[k]] ≤ halfQ ? Bred(a[1][idxa[k]], P[i]) : Bred(signed(a[1][idxa[k]] - Q[1].Q), P[i])
                end
            end
        end
    else
        Qtilde, Qstar, QmodP, oneoverQ, isPinQ, v128, v64, buff = be.Qtilde, be.Qstar, be.QmodP, be.oneoverQ, be.isPinQ, be.v128, be.v64, be.buff

        @inbounds for k = 0:8:N-8
            for i = 1:Qlen
                buff[1, i] = Bmul(a[i][idxa[k+1]], Qtilde[i], Q[i])
                buff[2, i] = Bmul(a[i][idxa[k+2]], Qtilde[i], Q[i])
                buff[3, i] = Bmul(a[i][idxa[k+3]], Qtilde[i], Q[i])
                buff[4, i] = Bmul(a[i][idxa[k+4]], Qtilde[i], Q[i])
                buff[5, i] = Bmul(a[i][idxa[k+5]], Qtilde[i], Q[i])
                buff[6, i] = Bmul(a[i][idxa[k+6]], Qtilde[i], Q[i])
                buff[7, i] = Bmul(a[i][idxa[k+7]], Qtilde[i], Q[i])
                buff[8, i] = Bmul(a[i][idxa[k+8]], Qtilde[i], Q[i])
            end
            v128[1] = 0
            v128[2] = 0
            v128[3] = 0
            v128[4] = 0
            v128[5] = 0
            v128[6] = 0
            v128[7] = 0
            v128[8] = 0

            for i = 1:Qlen
                v128[1] += mult_and_round(buff[1, i], oneoverQ[i])
                v128[2] += mult_and_round(buff[2, i], oneoverQ[i])
                v128[3] += mult_and_round(buff[3, i], oneoverQ[i])
                v128[4] += mult_and_round(buff[4, i], oneoverQ[i])
                v128[5] += mult_and_round(buff[5, i], oneoverQ[i])
                v128[6] += mult_and_round(buff[6, i], oneoverQ[i])
                v128[7] += mult_and_round(buff[7, i], oneoverQ[i])
                v128[8] += mult_and_round(buff[8, i], oneoverQ[i])
            end
            v64[1] = round_to_uint64(v128[1])
            v64[2] = round_to_uint64(v128[2])
            v64[3] = round_to_uint64(v128[3])
            v64[4] = round_to_uint64(v128[4])
            v64[5] = round_to_uint64(v128[5])
            v64[6] = round_to_uint64(v128[6])
            v64[7] = round_to_uint64(v128[7])
            v64[8] = round_to_uint64(v128[8])

            for j = 1:Plen
                # Decompose more efficiently.
                if isPinQ[j] > 0
                    res[j][idxres[k+1]] = a[isPinQ[j]][idxa[k+1]]
                    res[j][idxres[k+2]] = a[isPinQ[j]][idxa[k+2]]
                    res[j][idxres[k+3]] = a[isPinQ[j]][idxa[k+3]]
                    res[j][idxres[k+4]] = a[isPinQ[j]][idxa[k+4]]
                    res[j][idxres[k+5]] = a[isPinQ[j]][idxa[k+5]]
                    res[j][idxres[k+6]] = a[isPinQ[j]][idxa[k+6]]
                    res[j][idxres[k+7]] = a[isPinQ[j]][idxa[k+7]]
                    res[j][idxres[k+8]] = a[isPinQ[j]][idxa[k+8]]
                    continue
                end

                res[j][idxres[k+1]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[1], P[j])
                res[j][idxres[k+2]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[2], P[j])
                res[j][idxres[k+3]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[3], P[j])
                res[j][idxres[k+4]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[4], P[j])
                res[j][idxres[k+5]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[5], P[j])
                res[j][idxres[k+6]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[6], P[j])
                res[j][idxres[k+7]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[7], P[j])
                res[j][idxres[k+8]] = lazy_Bmul(neg(QmodP[j], P[j]), v64[8], P[j])
                for i = 1:Qlen-1
                    res[j][idxres[k+1]] = lazy_Bred(widemul(buff[1, i], Qstar[i, j]) + res[j][idxres[k+1]], P[j])
                    res[j][idxres[k+2]] = lazy_Bred(widemul(buff[2, i], Qstar[i, j]) + res[j][idxres[k+2]], P[j])
                    res[j][idxres[k+3]] = lazy_Bred(widemul(buff[3, i], Qstar[i, j]) + res[j][idxres[k+3]], P[j])
                    res[j][idxres[k+4]] = lazy_Bred(widemul(buff[4, i], Qstar[i, j]) + res[j][idxres[k+4]], P[j])
                    res[j][idxres[k+5]] = lazy_Bred(widemul(buff[5, i], Qstar[i, j]) + res[j][idxres[k+5]], P[j])
                    res[j][idxres[k+6]] = lazy_Bred(widemul(buff[6, i], Qstar[i, j]) + res[j][idxres[k+6]], P[j])
                    res[j][idxres[k+7]] = lazy_Bred(widemul(buff[7, i], Qstar[i, j]) + res[j][idxres[k+7]], P[j])
                    res[j][idxres[k+8]] = lazy_Bred(widemul(buff[8, i], Qstar[i, j]) + res[j][idxres[k+8]], P[j])
                end
                res[j][idxres[k+1]] = Bred(widemul(buff[1, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+1]], P[j])
                res[j][idxres[k+2]] = Bred(widemul(buff[2, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+2]], P[j])
                res[j][idxres[k+3]] = Bred(widemul(buff[3, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+3]], P[j])
                res[j][idxres[k+4]] = Bred(widemul(buff[4, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+4]], P[j])
                res[j][idxres[k+5]] = Bred(widemul(buff[5, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+5]], P[j])
                res[j][idxres[k+6]] = Bred(widemul(buff[6, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+6]], P[j])
                res[j][idxres[k+7]] = Bred(widemul(buff[7, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+7]], P[j])
                res[j][idxres[k+8]] = Bred(widemul(buff[8, Qlen], Qstar[Qlen, j]) + res[j][idxres[k+8]], P[j])
            end
        end

        @inbounds for k = 8(N>>3)+1:N
            for i = 1:Qlen
                buff[i] = Bmul(a[i][idxa[k]], Qtilde[i], Q[i])
            end

            v128 = zero(UInt128)
            for i = 1:Qlen
                v128 += mult_and_round(buff[i], oneoverQ[i])
            end
            v64 = round_to_uint64(v128)

            for j = 1:Plen
                # Decompose more efficiently.
                if isPinQ[j] > 0
                    res[j][idxres[k]] = a[isPinQ[j]][idxa[k]]
                    continue
                end

                res[j][idxres[k]] = lazy_Bmul(neg(QmodP[j], P[j]), v64, P[j])
                for i = 1:Qlen-1
                    res[j][idxres[k]] = lazy_Bred(widemul(buff[i], Qstar[i, j]) + res[j][idxres[k]], P[j])
                end
                res[j][idxres[k]] = Bred(widemul(buff[Qlen], Qstar[Qlen, j]) + res[j][idxres[k]], P[j])
            end
        end
    end

    return nothing
end

"""
SimpleScaler provides a simple scaling operation from P to Q.
To be precise, it computes x (mod P) → ⌊Q/P ⋅ x⌉ (mod Q).
"""
struct SimpleScaler
    P::Vector{Modulus}
    Q::Vector{Modulus}
    Ptilde::Vector{UInt64}
    ω::Array{UInt64,2}
    θ::Vector{UInt128}
    v128::Vector{UInt128}
    buffP::Array{UInt64,2}

    @views function SimpleScaler(P::Moduli, Q::Moduli)::SimpleScaler
        Plen, Qlen = length(P), length(Q)

        # Compute required values.
        Ptilde = Vector{UInt64}(undef, Plen)
        ω = Array{UInt64}(undef, Plen, Qlen)
        θ = Vector{UInt128}(undef, Plen)
        v128 = Vector{UInt128}(undef, 8)
        buffP = Array{UInt64}(undef, 8, Plen)

        # Compute Rtilde
        @inbounds for i = 1:Plen
            Ptilde[i] = 1
            for j = 1:Plen
                i == j && continue
                Ptilde[i] = Bmul(Ptilde[i], invmod(P[j].Q, P[i].Q), P[i])
            end
        end

        # Compute ω and θ
        setprecision(192)
        Qbig = prod(Q)
        @inbounds for i = 1:Plen
            ωi = Qbig ÷ P[i].Q
            for j = 1:Qlen
                ω[i, j] = (ωi % Q[j].Q) % UInt64
            end
            θ[i] = round(UInt128, (Int128(1) << fixed_prec) * (Qbig - ωi * P[i].Q) / P[i].Q)
        end

        new(collect(P), collect(Q), Ptilde, ω, θ, v128, buffP)
    end
end

@views function simple_scale_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, ss::SimpleScaler)::Nothing
    N = length(res[1])
    P, Q = ss.P, ss.Q
    Plen, Qlen = length(P), length(Q)

    if length(res[1]) ≠ length(a[1]) ≠ N
        throw(DimensionMismatch("The input and output arrays do not match."))
    end
    if length(res) ≠ Qlen && length(a) ≠ Plen
        throw(DimensionMismatch("The input and output arrays do not match."))
    end

    Ptilde, ω, θ = ss.Ptilde, ss.ω, ss.θ
    v128, buffP = ss.v128, ss.buffP

    v128[1] = zero(UInt128)
    @inbounds for i = 1:Plen
        buffP[1, i] = Bmul(a[i], Ptilde[i], P[i])
        v128[1] += mult_and_round(buffP[1, i], θ[i])
    end
    v128[1] = round_to_uint128(v128[1])

    @inbounds for i = 1:Qlen
        res[i] = lazy_Bred(v128[1], Q[i])
        for j = 1:Plen-1
            res[i] = lazy_Bred(widemul(buffP[1, j], ω[j, i]) + res[i], Q[i])
        end
        res[i] = Bred(widemul(buffP[1, end], ω[end, i]) + res[i], Q[i])
    end

    return nothing
end

@views function simple_scale_to!(res::AbstractVector{Vector{UInt64}}, a::AbstractVector{Vector{UInt64}}, ss::SimpleScaler)::Nothing
    N = length(res[1])
    P, Q = ss.P, ss.Q
    Plen, Qlen = length(P), length(Q)

    if length(res[1]) ≠ length(a[1]) ≠ N
        throw(DimensionMismatch("The input and output arrays do not match."))
    end
    if length(res) ≠ Qlen && length(a) ≠ Plen
        throw(DimensionMismatch("The input and output arrays do not match."))
    end

    Ptilde, ω, θ = ss.Ptilde, ss.ω, ss.θ
    v128, buffP = ss.v128, ss.buffP

    @inbounds for k = 0:8:N-8
        v128[1] = zero(UInt128)
        v128[2] = zero(UInt128)
        v128[3] = zero(UInt128)
        v128[4] = zero(UInt128)
        v128[5] = zero(UInt128)
        v128[6] = zero(UInt128)
        v128[7] = zero(UInt128)
        v128[8] = zero(UInt128)

        for i = 1:Plen
            buffP[1, i] = Bmul(a[i][k+1], Ptilde[i], P[i])
            buffP[2, i] = Bmul(a[i][k+2], Ptilde[i], P[i])
            buffP[3, i] = Bmul(a[i][k+3], Ptilde[i], P[i])
            buffP[4, i] = Bmul(a[i][k+4], Ptilde[i], P[i])
            buffP[5, i] = Bmul(a[i][k+5], Ptilde[i], P[i])
            buffP[6, i] = Bmul(a[i][k+6], Ptilde[i], P[i])
            buffP[7, i] = Bmul(a[i][k+7], Ptilde[i], P[i])
            buffP[8, i] = Bmul(a[i][k+8], Ptilde[i], P[i])
            v128[1] += mult_and_round(buffP[1, i], θ[i])
            v128[2] += mult_and_round(buffP[2, i], θ[i])
            v128[3] += mult_and_round(buffP[3, i], θ[i])
            v128[4] += mult_and_round(buffP[4, i], θ[i])
            v128[5] += mult_and_round(buffP[5, i], θ[i])
            v128[6] += mult_and_round(buffP[6, i], θ[i])
            v128[7] += mult_and_round(buffP[7, i], θ[i])
            v128[8] += mult_and_round(buffP[8, i], θ[i])
        end
        v128[1] = round_to_uint128(v128[1])
        v128[2] = round_to_uint128(v128[2])
        v128[3] = round_to_uint128(v128[3])
        v128[4] = round_to_uint128(v128[4])
        v128[5] = round_to_uint128(v128[5])
        v128[6] = round_to_uint128(v128[6])
        v128[7] = round_to_uint128(v128[7])
        v128[8] = round_to_uint128(v128[8])

        for i = 1:Qlen
            res[i][k+1] = lazy_Bred(v128[1], Q[i])
            res[i][k+2] = lazy_Bred(v128[2], Q[i])
            res[i][k+3] = lazy_Bred(v128[3], Q[i])
            res[i][k+4] = lazy_Bred(v128[4], Q[i])
            res[i][k+5] = lazy_Bred(v128[5], Q[i])
            res[i][k+6] = lazy_Bred(v128[6], Q[i])
            res[i][k+7] = lazy_Bred(v128[7], Q[i])
            res[i][k+8] = lazy_Bred(v128[8], Q[i])
            for j = 1:Plen-1
                res[i][k+1] = lazy_Bred(widemul(buffP[1, j], ω[j, i]) + res[i][k+1], Q[i])
                res[i][k+2] = lazy_Bred(widemul(buffP[2, j], ω[j, i]) + res[i][k+2], Q[i])
                res[i][k+3] = lazy_Bred(widemul(buffP[3, j], ω[j, i]) + res[i][k+3], Q[i])
                res[i][k+4] = lazy_Bred(widemul(buffP[4, j], ω[j, i]) + res[i][k+4], Q[i])
                res[i][k+5] = lazy_Bred(widemul(buffP[5, j], ω[j, i]) + res[i][k+5], Q[i])
                res[i][k+6] = lazy_Bred(widemul(buffP[6, j], ω[j, i]) + res[i][k+6], Q[i])
                res[i][k+7] = lazy_Bred(widemul(buffP[7, j], ω[j, i]) + res[i][k+7], Q[i])
                res[i][k+8] = lazy_Bred(widemul(buffP[8, j], ω[j, i]) + res[i][k+8], Q[i])
            end
            res[i][k+1] = Bred(widemul(buffP[1, end], ω[end, i]) + res[i][k+1], Q[i])
            res[i][k+2] = Bred(widemul(buffP[2, end], ω[end, i]) + res[i][k+2], Q[i])
            res[i][k+3] = Bred(widemul(buffP[3, end], ω[end, i]) + res[i][k+3], Q[i])
            res[i][k+4] = Bred(widemul(buffP[4, end], ω[end, i]) + res[i][k+4], Q[i])
            res[i][k+5] = Bred(widemul(buffP[5, end], ω[end, i]) + res[i][k+5], Q[i])
            res[i][k+6] = Bred(widemul(buffP[6, end], ω[end, i]) + res[i][k+6], Q[i])
            res[i][k+7] = Bred(widemul(buffP[7, end], ω[end, i]) + res[i][k+7], Q[i])
            res[i][k+8] = Bred(widemul(buffP[8, end], ω[end, i]) + res[i][k+8], Q[i])
        end
    end

    @inbounds for k = 8(N>>3)+1:N
        v128[1] = zero(UInt128)
        for i = 1:Plen
            buffP[1, i] = Bmul(a[i][k], Ptilde[i], P[i])
            v128[1] += mult_and_round(buffP[1, i], θ[i])
        end
        v128[1] = round_to_uint128(v128[1])

        for i = 1:Qlen
            res[i][k] = lazy_Bred(v128[1], Q[i])
            for j = 1:Plen-1
                res[i][k] = lazy_Bred(widemul(buffP[1, j], ω[j, i]) + res[i][k], Q[i])
            end
            res[i][k] = Bred(widemul(buffP[1, end], ω[end, i]) + res[i][k], Q[i])
        end
    end

    return nothing
end

"""
    ComplexScaler provides a complex scaling operation from P to Q using an optimised HPS18-like algorithm.
    More precisely, it computes x (mod P) -> ⌊x ⋅ scale⌉ (mod Q).
"""
struct ComplexScaler
    P::Vector{Modulus}
    Q::Vector{Modulus}
    Ptilde::Vector{UInt64}
    oneoverP::Vector{UInt128}
    ω::Array{UInt64,2}
    θ::Vector{UInt128}
    ζ::Vector{UInt64}
    λ::UInt128
    v128::Vector{UInt128}
    v64::Vector{UInt64}
    buffP::Array{UInt64,2}

    function ComplexScaler(P::Moduli, Q::Moduli, scale::Real)::ComplexScaler
        Plen, Qlen = length(P), length(Q)

        Ptilde = Vector{UInt64}(undef, Plen)
        oneoverP = Vector{UInt128}(undef, Plen)
        ω = Array{UInt64}(undef, Plen, Qlen)
        θ = Vector{UInt128}(undef, Plen)
        ζ = Vector{UInt64}(undef, Qlen)
        v128 = Vector{UInt128}(undef, 8)
        v64 = Vector{UInt64}(undef, 8)
        buffP = Array{UInt64}(undef, 8, Plen)

        Pbig = prod(P)

        # Compute Ptilde, ωi, θi, oneoverP.
        @inbounds for i = 1:Plen
            Ptilde[i] = 1
            for j = 1:Plen
                i == j && continue
                Ptilde[i] = Bmul(Ptilde[i], invmod(P[j].Q, P[i].Q), P[i])
            end

            Ptildei = Pbig ÷ P[i].Q
            tmp = Ptildei * scale
            ωi = floor(BigInt, tmp)
            for j = 1:Qlen
                ω[i, j] = (ωi % Q[j].Q) % UInt64
            end
            θ[i] = round(UInt128, (Int128(1) << fixed_prec) * (tmp - ωi))
            oneoverP[i] = round(UInt128, (Int128(1) << fixed_prec) / big(P[i].Q))
        end

        # Compute ζ and λ.
        tmp = scale * Pbig
        ζbig = floor(BigInt, tmp)
        λ = round(UInt128, (Int128(1) << fixed_prec) * (tmp - ζbig))
        for i = 1:Qlen
            ζ[i] = ζbig % Q[i].Q
        end

        new(collect(P), collect(Q), Ptilde, oneoverP, ω, θ, ζ, λ, v128, v64, buffP)
    end
end

@views function complex_scale_to!(res::AbstractVector{Vector{UInt64}}, a::AbstractVector{Vector{UInt64}}, cs::ComplexScaler)::Nothing
    N = length(res[1])
    P, Q = cs.P, cs.Q
    Plen, Qlen = length(P), length(Q)

    if length(a[1]) ≠ length(res[1]) ≠ N
        throw(DimensionMismatch("The input and output arrays do not match."))
    end
    if length(res) ≠ Qlen && length(a) ≠ Plen
        throw(DimensionMismatch("The input and output arrays do not match."))
    end

    Ptilde, oneoverP, ω, θ, ζ, λ = cs.Ptilde, cs.oneoverP, cs.ω, cs.θ, cs.ζ, cs.λ
    v128, v64, buffP = cs.v128, cs.v64, cs.buffP

    @inbounds for k = 0:8:N-8
        v128[1] = zero(UInt128)
        v128[2] = zero(UInt128)
        v128[3] = zero(UInt128)
        v128[4] = zero(UInt128)
        v128[5] = zero(UInt128)
        v128[6] = zero(UInt128)
        v128[7] = zero(UInt128)
        v128[8] = zero(UInt128)

        for i = 1:Plen
            buffP[1, i] = Bmul(a[i][k+1], Ptilde[i], P[i])
            buffP[2, i] = Bmul(a[i][k+2], Ptilde[i], P[i])
            buffP[3, i] = Bmul(a[i][k+3], Ptilde[i], P[i])
            buffP[4, i] = Bmul(a[i][k+4], Ptilde[i], P[i])
            buffP[5, i] = Bmul(a[i][k+5], Ptilde[i], P[i])
            buffP[6, i] = Bmul(a[i][k+6], Ptilde[i], P[i])
            buffP[7, i] = Bmul(a[i][k+7], Ptilde[i], P[i])
            buffP[8, i] = Bmul(a[i][k+8], Ptilde[i], P[i])
            v128[1] += mult_and_round(buffP[1, i], oneoverP[i])
            v128[2] += mult_and_round(buffP[2, i], oneoverP[i])
            v128[3] += mult_and_round(buffP[3, i], oneoverP[i])
            v128[4] += mult_and_round(buffP[4, i], oneoverP[i])
            v128[5] += mult_and_round(buffP[5, i], oneoverP[i])
            v128[6] += mult_and_round(buffP[6, i], oneoverP[i])
            v128[7] += mult_and_round(buffP[7, i], oneoverP[i])
            v128[8] += mult_and_round(buffP[8, i], oneoverP[i])
        end
        v64[1] = round_to_uint64(v128[1])
        v64[2] = round_to_uint64(v128[2])
        v64[3] = round_to_uint64(v128[3])
        v64[4] = round_to_uint64(v128[4])
        v64[5] = round_to_uint64(v128[5])
        v64[6] = round_to_uint64(v128[6])
        v64[7] = round_to_uint64(v128[7])
        v64[8] = round_to_uint64(v128[8])

        for i = 1:Qlen
            res[i][k+1] = neg(Bmul(v64[1], ζ[i], Q[i]), Q[i])
            res[i][k+2] = neg(Bmul(v64[2], ζ[i], Q[i]), Q[i])
            res[i][k+3] = neg(Bmul(v64[3], ζ[i], Q[i]), Q[i])
            res[i][k+4] = neg(Bmul(v64[4], ζ[i], Q[i]), Q[i])
            res[i][k+5] = neg(Bmul(v64[5], ζ[i], Q[i]), Q[i])
            res[i][k+6] = neg(Bmul(v64[6], ζ[i], Q[i]), Q[i])
            res[i][k+7] = neg(Bmul(v64[7], ζ[i], Q[i]), Q[i])
            res[i][k+8] = neg(Bmul(v64[8], ζ[i], Q[i]), Q[i])
            for j = 1:Plen
                res[i][k+1] = Bred(widemul(buffP[1, j], ω[j, i]) + res[i][k+1], Q[i])
                res[i][k+2] = Bred(widemul(buffP[2, j], ω[j, i]) + res[i][k+2], Q[i])
                res[i][k+3] = Bred(widemul(buffP[3, j], ω[j, i]) + res[i][k+3], Q[i])
                res[i][k+4] = Bred(widemul(buffP[4, j], ω[j, i]) + res[i][k+4], Q[i])
                res[i][k+5] = Bred(widemul(buffP[5, j], ω[j, i]) + res[i][k+5], Q[i])
                res[i][k+6] = Bred(widemul(buffP[6, j], ω[j, i]) + res[i][k+6], Q[i])
                res[i][k+7] = Bred(widemul(buffP[7, j], ω[j, i]) + res[i][k+7], Q[i])
                res[i][k+8] = Bred(widemul(buffP[8, j], ω[j, i]) + res[i][k+8], Q[i])
            end
        end

        v128[1] = zero(UInt128)
        v128[2] = zero(UInt128)
        v128[3] = zero(UInt128)
        v128[4] = zero(UInt128)
        v128[5] = zero(UInt128)
        v128[6] = zero(UInt128)
        v128[7] = zero(UInt128)
        v128[8] = zero(UInt128)
        for i = 1:Plen
            v128[1] += mult_and_round(buffP[1, i], θ[i])
            v128[2] += mult_and_round(buffP[2, i], θ[i])
            v128[3] += mult_and_round(buffP[3, i], θ[i])
            v128[4] += mult_and_round(buffP[4, i], θ[i])
            v128[5] += mult_and_round(buffP[5, i], θ[i])
            v128[6] += mult_and_round(buffP[6, i], θ[i])
            v128[7] += mult_and_round(buffP[7, i], θ[i])
            v128[8] += mult_and_round(buffP[8, i], θ[i])
        end
        v128[1] -= mult_and_round(v64[1], λ)
        v128[2] -= mult_and_round(v64[2], λ)
        v128[3] -= mult_and_round(v64[3], λ)
        v128[4] -= mult_and_round(v64[4], λ)
        v128[5] -= mult_and_round(v64[5], λ)
        v128[6] -= mult_and_round(v64[6], λ)
        v128[7] -= mult_and_round(v64[7], λ)
        v128[8] -= mult_and_round(v64[8], λ)

        v128[1] = round_to_uint128(signed(v128[1]))
        v128[2] = round_to_uint128(signed(v128[2]))
        v128[3] = round_to_uint128(signed(v128[3]))
        v128[4] = round_to_uint128(signed(v128[4]))
        v128[5] = round_to_uint128(signed(v128[5]))
        v128[6] = round_to_uint128(signed(v128[6]))
        v128[7] = round_to_uint128(signed(v128[7]))
        v128[8] = round_to_uint128(signed(v128[8]))

        for i = 1:Qlen
            res[i][k+1] = Bred(signed(v128[1]) + res[i][k+1], Q[i])
            res[i][k+2] = Bred(signed(v128[2]) + res[i][k+2], Q[i])
            res[i][k+3] = Bred(signed(v128[3]) + res[i][k+3], Q[i])
            res[i][k+4] = Bred(signed(v128[4]) + res[i][k+4], Q[i])
            res[i][k+5] = Bred(signed(v128[5]) + res[i][k+5], Q[i])
            res[i][k+6] = Bred(signed(v128[6]) + res[i][k+6], Q[i])
            res[i][k+7] = Bred(signed(v128[7]) + res[i][k+7], Q[i])
            res[i][k+8] = Bred(signed(v128[8]) + res[i][k+8], Q[i])
        end
    end

    @inbounds for k = 8(N>>3)+1:N
        v128[1] = zero(UInt128)
        for i = 1:Plen
            buffP[1, i] = Bmul(a[i][k], Ptilde[i], P[i])
            v128[1] += mult_and_round(buffP[1, i], oneoverP[i])
        end
        v64[1] = round_to_uint64(v128[1])

        for i = 1:Qlen
            res[i][k] = neg(Bmul(v64[1], ζ[i], Q[i]), Q[i])
            for j = 1:Plen
                res[i][k] = Bred(widemul(buffP[1, j], ω[j, i]) + res[i][k], Q[i])
            end
        end

        v128[1] = zero(UInt128)
        for i = 1:Plen
            v128[1] += mult_and_round(buffP[1, i], θ[i])
        end
        v128[1] -= mult_and_round(v64[1], λ)
        v128[1] = round_to_uint128(signed(v128[1]))

        for i = 1:Qlen
            res[i][k] = Bred(signed(v128[1]) + res[i][k], Q[i])
        end
    end

    return nothing
end