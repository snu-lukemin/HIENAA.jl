"""
This function evaluates the automorphism X -> Xⁱᵈˣ.
"""
@views function _automorphism!(idx::Int64, a::AbstractVector{UInt64}, isntt::Bool, autbuff::AbstractVector{UInt64}, Q::Modulus)::Nothing
    N = length(autbuff)
    m = N << 1

    if idx ≥ 0
        idx = idx & (m - 1)
    else
        idx = m - ((-idx) & (m - 1))
    end

    if isntt
        @. autbuff = a
        scramble!(autbuff, 2)
        @inbounds @simd for i = 1:N
            a[i] = autbuff[(((2i-1)*idx)&(m-1)+1)>>1]
        end
        scramble!(a, 2)
    else
        @. autbuff = a
        @inbounds for i = 0:N-1
            j = (idx * i) & (m - 1)
            if j < N
                a[j+1] = autbuff[i+1]
            else
                a[j-N+1] = neg(autbuff[i+1], Q)
            end
        end
    end

    return nothing
end

"""
This function evaluates the automorphism X -> Xⁱᵈˣ.
"""
@views function _automorphism!(idx::Int64, a::AbstractVector{UInt64}, isntt::Bool, autbuff::Vector{UInt64}, ntter::CyclotomicNTTransformerArb, rdtor::ReductorCyclo)::Nothing
    m, N = length(autbuff), length(a)
    idx = Base.mod(idx, m)

    if isntt
        autidx = ntter.autidxset[idx+1]
        dims = ntter.dims

        if length(dims) == 1
            j1 = autidx[1]
            j2 = 0
            @inbounds for _ = 1:dims[1]
                autbuff[j2+1] = a[j1%N+1]
                j1 += 1
                j2 += 1
            end
        elseif length(dims) == 2
            j1 = autidx[1] + autidx[2] * dims[1]
            j2 = 0
            @inbounds for i1 = 1:dims[1]
                for _ = 1:dims[2]
                    autbuff[j2+1] = a[j1%N+1]
                    j1 += dims[1]
                    j2 += dims[1]
                end
                i1 + autidx[1] == dims[1] && (j1 -= dims[1])
                j1 += 1
                j2 += 1 - N
            end
        elseif length(dims) == 3
            dim12 = dims[1] * dims[2]

            j1 = autidx[1] + autidx[2] * dims[1] + autidx[3] * dim12
            j2 = 0
            @inbounds for i1 = 1:dims[1]
                for i2 = 1:dims[2]
                    for _ = 1:dims[3]
                        autbuff[j2+1] = a[j1%N+1]
                        j1 += dim12
                        j2 += dim12
                    end
                    i2 + autidx[2] == dims[2] && (j1 -= dim12)
                    j1 += dims[1] - N
                    j2 += dims[1] - N
                end
                i1 + autidx[1] == dims[1] && (j1 -= dims[1])
                j1 += 1
                j2 += 1 - dim12
            end
        else
            dim12, dim123 = dims[1] * dims[2], dims[1] * dims[2] * dims[3]

            j1 = autidx[1] + autidx[2] * dims[1] + autidx[3] * dim12 + autidx[4] * dim123
            j2 = 0
            @inbounds for i1 = 1:dims[1]
                for i2 = 1:dims[2]
                    for i3 = 1:dims[3]
                        for _ = 1:dims[4]
                            autbuff[j2+1] = a[j1%N+1]
                            j1 += dim123
                            j2 += dim123
                        end
                        i3 + autidx[3] == dims[3] && (j1 -= dim123)
                        j1 += dim12 - N
                        j2 += dim12 - N
                    end
                    i2 + autidx[2] == dims[2] && (j1 -= dim12)
                    j1 += dims[1]
                    j2 += dims[1] - dim123
                end
                i1 + autidx[1] == dims[1] && (j1 -= dims[1])
                j1 += 1
                j2 += 1 - dim12
            end
        end

        @inbounds for i = 1:N
            a[i] = autbuff[i]
        end
    else
        @. autbuff = 0
        autbuff[1] = a[1]
        @inbounds @simd for i = 1:N-1
            autbuff[(idx*i)%m+1] = a[i+1]
        end

        reduce!(autbuff, rdtor)

        @. a = autbuff[1:N]
    end

    return nothing
end