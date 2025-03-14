function get_c2s_matrix(oper::BGVOperator; packlen::Int64=0)::PlainMatrix
    if ismissing(oper.packer)
        throw(ErrorException("The operator must have a packer."))
    end
    if !isa(oper.packer, IntPackerSubring)
        throw(ErrorException("The packer must be a subring packer."))
    end

    packlen == 0 && (packlen = oper.packer.k)
    if oper.packer.k % packlen ≠ 0
        throw(DomainError("The number of slots must divide the packing parameter."))
    end

    m, p, r, pr = oper.packer.m, oper.packer.p, oper.packer.r, oper.packer.pr

    resol = load_resolution(m, p, r)
    Bred_to!(resol, resol, pr)
    resollen = length(resol)

    if packlen < resollen
        @inbounds for i = 1:resollen÷packlen-1
            @views add_to!(resol[1:packlen], resol[1:packlen], resol[i*packlen+1:(i+1)*packlen], oper.packer.pr)
        end
        resize!(resol, packlen)

        resollen = packlen
    end

    level = length(oper.Qatlevel) - 1
    res = Dict{Int64,PlainPoly}()
    mattmp = Vector{UInt64}(undef, packlen)
    n1, n2 = get_bsgs_param(packlen)

    @inbounds for i = 0:n2-1
        for j = 0:n1-1
            i * n1 + j ≥ packlen && break

            @inbounds for k = 0:packlen-1
                idx = mod(-2k + i * n1 + j, packlen) + 1
                mattmp[k+1] = resol[(packlen-idx+1)%packlen+1]
            end
            all(iszero, mattmp) && continue
            circshift!(mattmp, -i * n1)
            res[i*n1+j] = encode(mattmp, oper, level=level)
        end
    end

    PlainMatrix(res, (n1, n2))
end

function get_s2c_matrix(oper::BGVOperator; packlen::Int64=0)::PlainMatrix
    if ismissing(oper.packer)
        throw(ErrorException("The operator must have a packer."))
    end
    if !isa(oper.packer, IntPackerSubring)
        throw(ErrorException("The packer must be a subring packer."))
    end

    packlen == 0 && (packlen = oper.packer.k)
    if oper.packer.k % packlen ≠ 0
        throw(DomainError("The number of slots must divide the packing parameter."))
    end

    m, p, r, pr = oper.packer.m, oper.packer.p, oper.packer.r, oper.packer.pr

    resol = load_resolution(m, p, r)
    Bred_to!(resol, resol, pr)
    resollen = length(resol)

    if packlen < resollen
        @inbounds for i = 1:resollen÷packlen-1
            @views add_to!(resol[1:packlen], resol[1:packlen], resol[i*packlen+1:(i+1)*packlen], oper.packer.pr)
        end
        resize!(resol, packlen)

        resollen = packlen
    end
    ordp = (m - 1) ÷ resollen

    isodd(ordp) && circshift!(resol, resollen >> 1)
    Bmul_to!(resol, UInt64(m), resol, pr)
    add_to!(resol, resol, Bred(ordp, pr), pr)

    level = length(oper.Qatlevel) - 1
    res = Dict{Int64,PlainPoly}()
    mattmp = Vector{UInt64}(undef, packlen)
    n1, n2 = get_bsgs_param(packlen)

    @inbounds for i = 0:n2-1
        for j = 0:n1-1
            i * n1 + j ≥ packlen && break

            for k = 0:packlen-1
                idx = mod(-2k + i * n1 + j, packlen) + 1
                mattmp[k+1] = resol[(packlen-idx+1)%packlen+1]
            end
            all(iszero, mattmp) && continue
            circshift!(mattmp, -i * n1)
            res[i*n1+j] = encode(mattmp, oper, level=level)
        end
    end

    PlainMatrix(res, (n1, n2))
end

function get_digit_extraction_poly(i::Int64, interp::Interpolator)::Vector{UInt64}
    p, r = interp.p, interp.r

    if p == 2
        x = UInt64.(collect(0:i+1))
        y = @. x % p

        res = interpolate(x, y, interp)
        @views @. res[2:2:end] = 0
    else
        x = UInt64.(collect(0:(p-1)*(i-1)+1))
        y = @. mod((Int64(x) + p ÷ 2) % p - p ÷ 2, p^r)

        res = interpolate(x, y, interp)
        @views @. res[1:2:end] = 0
    end

    res
end