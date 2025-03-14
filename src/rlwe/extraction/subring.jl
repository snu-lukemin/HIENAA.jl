"""
An Extractor type for the subring.
"""
struct ExtractorSubring <: Extractor
    m::Int64
    N::Int64
    d::Int64
    idx::Vector{Int64}

    function ExtractorSubring(param::SubringParam)::ExtractorSubring
        m, N = param.m, param.N
        d, g = param.d, primitive_root_finder(m)

        idx = Vector{Int64}(undef, m - 1)
        @inbounds for i = 0:m-2
            idx[powermod(g, i, m)] = i % N + 1
        end

        new(m, N, d, idx)
    end
end

"""
extract performs sample extraction for input RLWE ciphertext, and stores the value in res.
"""
function extract(ct::RLWE, Q::Modulus, extor::ExtractorSubring; rowidx::Union{AbstractRange{Int64},Missing}=missing, colidx::Union{AbstractRange{Int64},Missing}=missing)::Tuple{Matrix{UInt64}, Vector{UInt64}}
    if length(ct) ≠ 1
        throw(DomainError("The input ciphertext must have a single level."))
    end
    if ct.b.isntt[] || ct.a.isntt[]
        throw(DomainError("The input ciphertext must be in coefficient domain."))
    end

    m, N, idx = extor.m, extor.N, extor.idx

    ismissing(rowidx) && (rowidx = 1:N)
    ismissing(colidx) && (colidx = 1:N)

    reslen, collen = length(rowidx), length(colidx)
    A, b = zeros(UInt64, reslen, collen), zeros(UInt64, reslen)

    g = primitive_root_finder(m)
    @inbounds for i = eachindex(rowidx)
        gi = powermod(g, rowidx[i] - 1, m)
        for j = 1:m-1
            jidx = findfirst(isequal(idx[j]), colidx)
            isnothing(jidx) && continue

            if gi == j
                tmp = neg(ct.a.coeffs[1][idx[m-j]], Q)
            elseif gi < j
                tmp = sub(ct.a.coeffs[1][idx[gi+m-j]], ct.a.coeffs[1][idx[m-j]], Q)
            else
                tmp = sub(ct.a.coeffs[1][idx[gi-j]], ct.a.coeffs[1][idx[m-j]], Q)
            end
           
            A[i, jidx] = add(A[i, jidx], tmp, Q)
        end
        b[i] = ct.b.coeffs[1][rowidx[i]]
    end

    A, b
end