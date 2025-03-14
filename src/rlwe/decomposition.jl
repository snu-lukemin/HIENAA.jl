"""
Decomposer is a struct for the gadget decomposition. It decomposes the polynomial using RNS decomposition, from bottom to top, with digit length `dlen`.
"""
struct Decomposer
    be::Vector{BasisExtender}
    gvec::Vector{ModScalar}
    Plen::Int64
    Qlen::Int64
    dlen::Int64
    glen::Int64

    function Decomposer(P::Moduli, Q::Moduli, dlen::Int64=0)::Decomposer
        Plen, Qlen = length(P), length(Q)
        dlen == 0 && (dlen = Plen)
        glen = ceil(Int64, Qlen / dlen)

        PQ = vcat(P, Q)
        bigP = prod(P)

        be = Vector{BasisExtender}(undef, glen)
        gvec = Vector{ModScalar}(undef, glen)
        @views @inbounds for i = 1:glen
            lo, hi = 1 + (i - 1) * dlen, min(i * dlen, Qlen)
            be[i] = BasisExtender(Q[lo:hi], PQ)

            tildeQ, Di = big(1), big(1)
            for j = 1:Qlen
                if lo ≤ j ≤ hi
                    Di *= Q[j].Q
                else
                    tildeQ *= Q[j].Q
                end
            end
            tildeQ *= invmod(tildeQ, Di) * bigP

            gveci = zeros(UInt64, Qlen + Plen)
            for j = 1:Qlen
                gveci[Plen+j] = (tildeQ % Q[j].Q) % UInt64
            end
            gvec[i] = ModScalar(gveci)
        end

        new(be, gvec, Plen, Qlen, dlen, glen)
    end

    function Decomposer(Q::Moduli, dlen::Int64=1)::Decomposer
        Qlen = length(Q)
        glen = ceil(Int64, Qlen / dlen)

        be = Vector{BasisExtender}(undef, glen)
        gvec = Vector{ModScalar}(undef, glen)
        @views @inbounds for i = 1:glen
            lo, hi = 1 + (i - 1) * dlen, min(i * dlen, Qlen)
            be[i] = BasisExtender(Q[lo:hi], Q)

            tildeQ, Di = big(1), big(1)
            for j = 1:Qlen
                if lo ≤ j ≤ hi
                    Di *= Q[j].Q
                else
                    tildeQ *= Q[j].Q
                end
            end
            tildeQ *= invmod(tildeQ, Di)

            gveci = zeros(UInt64, Qlen)
            for j = 1:Qlen
                gveci[j] = (tildeQ % Q[j].Q) % UInt64
            end
            gvec[i] = ModScalar(gveci)
        end

        new(be, gvec, 0, Qlen, dlen, glen)
    end
end

function decompose(x::PlainPoly, decer::Decomposer)::Tensor
    xlen, Plen, dlen = length(x), decer.Plen, decer.dlen
    res = Tensor(x.val.N, Plen + xlen, ceil(Int64, xlen / dlen))
    decompose_to!(res, x, decer)
    res
end

@views function decompose_to!(res::Tensor, x::PlainPoly, decer::Decomposer)::Nothing
    be, Plen, dlen, glen, xlen = decer.be, decer.Plen, decer.dlen, decer.glen, length(x)

    if x.val.isntt[]
        throw(DomainError("The input polynomial should not be in the evaluation form."))
    end
    if dlen * glen < xlen
        throw(DomainError("The decomposition parameter does not match the input polynomial."))
    end

    degree, len = size(res)
    if degree ≠ ceil(Int64, xlen / dlen) || len ≠ xlen + Plen
        throw(DomainError("The decomposition parameter does not match the input and output polynomials."))
    end

    isPQ, auxQ = (Plen ≠ 0), x.auxQ[]
    @inbounds for i = 1:degree
        idx1, idx2 = 1 + (i - 1) * dlen, min(i * dlen, xlen)
        basis_extend_to!(res.vals[i].coeffs, x.val.coeffs[idx1:idx2], be[i])
        res.vals[i].isntt[] = false
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end