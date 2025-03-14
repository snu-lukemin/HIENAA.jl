struct CKKSBootKey <: HEBootKey
    param::CKKSBootParameters
    ksktosparse::Union{RLEV,Nothing}
    kskfromsparse::Union{RLEV,Nothing}
    atk::Dict{Int64,RLEV}
    rlk::RLEV
end

function bootstrap_keygen(bootparam::CKKSBootParameters, oper::CKKSOperator, entor::SKEncryptor)::CKKSBootKey
    # Compute P0 and Q0.
    P0 = find_prime(oper.operQ.param, bootparam.P0_bits)[1]
    Q0 = find_prime(oper.operQ.param, bootparam.Q0_bits)[1]
    param = RLWEParameters(oper.operQ.param, UInt64[P0], UInt64[Q0], 1)
    operQ0 = Operator(param)

    # Generate sparse key.
    us = UniformSampler()
    sparse_hw = bootparam.sparse_hw

    if sparse_hw > 0
        sparsekey = ternary_ringkey(us, oper.operQ.param.N, sparse_hw)
        sparseenc = SKEncryptor(sparsekey, 3.2, operQ0)

        # Generate encapsulation keys.
        evalQ0 = operQ0.evalQ
        evalQL = oper.operQ.evalQ

        ksktosparse = rlev_encrypt(PlainPoly(RLWEkeyPQ(entor.key, evalQ0)), sparseenc)
        kskfromsparse = rlev_encrypt(PlainPoly(RLWEkeyPQ(sparsekey, evalQL)), entor)
    else
        ksktosparse = nothing
        kskfromsparse = nothing
    end

    # Set rotate keys
    c2s_radix, s2c_radix, packlen, N = bootparam.c2s_radix, bootparam.s2c_radix, bootparam.packlen, oper.operQ.param.N
    idxset = union(get_required_key_list_s2c(s2c_radix, packlen), get_required_key_list_c2s(c2s_radix, packlen))
    packer = oper.packer
    cube, cubegen = packer.cube, packer.cubegen
    atk = Dict{Int64,RLEV}()
    for i = eachindex(idxset)
        idx = idxset[i]
        autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
        atk[autidx] = automorphism_keygen(autidx, entor)
    end

    for i = trailing_zeros(packlen)+2:trailing_zeros(N)
        idx = (1 << i) + 1
        atk[idx] = automorphism_keygen(idx, entor)
    end

    # Generate conjugation key
    atk[-1] = automorphism_keygen(-1, entor)

    # Generate relinearisation key
    rlk = relin_keygen(entor)

    CKKSBootKey(bootparam, ksktosparse, kskfromsparse, atk, rlk)
end

struct CKKSBootstrapper <: HEBootstrapper
    param::CKKSBootParameters
    ksktosparse::Union{RLEV,Nothing}
    kskfromsparse::Union{RLEV,Nothing}
    C2Smatrix::Vector{PlainMatrix}
    S2Cmatrix::Vector{PlainMatrix}
    sinepoly::PolyCoeffs
    invpoly::PolyCoeffs
    operQ0::Operator
    operc2s::CKKSOperator
    opers2c::CKKSOperator
    opermod::CKKSOperator
    isthin::Bool

    function CKKSBootstrapper(param::CKKSBootParameters, ksktosparse::Union{RLEV,Nothing}, kskfromsparse::Union{RLEV,Nothing}, c2s::Vector{PlainMatrix}, s2c::Vector{PlainMatrix},
        sinepoly::PolyCoeffs, invpoly::PolyCoeffs, operQ0::Operator, operc2s::CKKSOperator, opers2c::CKKSOperator, opermod::CKKSOperator, isthin::Bool)::CKKSBootstrapper

        new(param, ksktosparse, kskfromsparse, c2s, s2c, sinepoly, invpoly, operQ0, operc2s, opers2c, opermod, isthin)
    end

    CKKSBootstrapper(bootkey::CKKSBootKey, oper::CKKSOperator)::CKKSBootstrapper = HEBootstrapper(bootkey, oper)
end

(::Type{HEBootstrapper})(bootkey::CKKSBootKey, oper::CKKSOperator)::CKKSBootstrapper = begin
    param, kskfromsparse, ksktosparse = bootkey.param, bootkey.kskfromsparse, bootkey.ksktosparse

    gap, packlen, Q0_bits = param.gap, param.packlen, param.Q0_bits
    c2s_radix, c2s_bits = param.c2s_radix, param.c2s_bits
    mod_bits, K, sine_fold, sine_degree, inverse_degree = param.mod_bits, param.K, param.sine_fold, param.sine_degree, param.inverse_degree
    s2c_radix, s2c_bits = param.s2c_radix, param.s2c_bits
    isthin = param.isthin

    # Generate operators
    P0 = find_prime(oper.operQ.param, param.P0_bits)[1]
    Q0 = find_prime(oper.operQ.param, param.Q0_bits)[1]
    operQ0 = Operator(RLWEParameters(oper.operQ.param, UInt64[P0], UInt64[Q0], 1))

    Qlen, auxQ = oper.Qatlevel[end]
    operc2s = CKKSOperator(oper, 2^c2s_bits, Qlen, auxQ)

    Qlen, auxQ = operc2s.Qatlevel[end-c2s_radix]
    opermod = CKKSOperator(oper, 2^mod_bits, Qlen, auxQ)

    if isthin
        # Find the top level for S2C.
        tarQ = Q0_bits + s2c_radix * s2c_bits

        Q, Qlen = oper.operQ.evalQ.moduli, 1
        local auxQ::UInt64
        while true
            nowQ = 0.0
            @inbounds for i = 1:Qlen
                nowQ += log2(Q[i].Q)
            end
            if nowQ ≥ tarQ
                auxQ = round(UInt64, 2^(tarQ - nowQ + log2(Q[Qlen].Q)))
                if auxQ == Q[Qlen].Q
                    auxQ = zero(UInt64)
                end
                break
            else
                Qlen += 1
            end
        end

        # Generate operators
        opers2c = CKKSOperator(oper, 2^s2c_bits, Qlen, auxQ)
    else
        evalmodlevel = trailing_zeros(Base._nextpow2(sine_degree + 1)) + trailing_zeros(Base._nextpow2(inverse_degree + 1))
        Qlen, auxQ = opermod.Qatlevel[end-evalmodlevel]
        opers2c = CKKSOperator(oper, 2^s2c_bits, Qlen, auxQ)
    end

    # Generate polynomial.
    m = 1 / gap
    offset = isodd(sine_fold) ? 0 : 1 / 4
    a, b = -K + 0.5 - offset, K - 0.5 + offset
    interval = [(BigFloat(i - m, precision=192), BigFloat(i + m, precision=192)) for i = -K+1:K-1]
    twopi = 2 * big(π)

    # Approximate the sine polynomial.
    if 2K - 1 > sine_degree
        throw(DomainError("The degree of the sine polynomial is too small."))
    end
    nc = approximate(x -> b * sinpi(2(x + offset) / sine_fold), interval, sine_degree, a, b)
    sinepoly = PolyCoeffs(nc, opermod, type=:chebyshev)

    # Compute the input interval of the inverse polynomial.
    interval = Tuple{BigFloat,BigFloat}[]
    for i ∈ 1:sine_fold
        lo, hi = sinpi(2big(i - m + offset) / sine_fold), sinpi(2big(i + m + offset) / sine_fold)
        if lo > hi
            lo, hi = hi, lo
        end
        push!(interval, (lo, hi))
    end
    sort!(interval, by=x -> x[1])

    # Approximate the inverse polynomial.
    if sine_fold > inverse_degree
        throw(DomainError("The degree of the inverse polynomial is too small."))
    end
    invnc = approximate(x -> asin(sin(sine_fold * asin(x) - twopi * offset)) / twopi, interval, inverse_degree, -1, 1)
    invpoly = PolyCoeffs(invnc, opermod, type=:chebyshev)

    # Generate matrices
    c2sscale = 0.5 * 2^mod_bits / (BigFloat(b) * Q0)
    s2cscale = (packlen < oper.operQ.param.N / 2 && isthin) ? 2.0 : 1.0

    c2s = get_c2s_matrix(operc2s, scale=c2sscale, radix=c2s_radix, packlen=packlen)
    s2c = get_s2c_matrix(opers2c, scale=s2cscale, radix=s2c_radix, packlen=packlen)

    CKKSBootstrapper(param, ksktosparse, kskfromsparse, c2s, s2c, sinepoly, invpoly, operQ0, operc2s, opers2c, opermod, isthin)
end

function modulus_after_bootstrap(boot::CKKSBootstrapper)::Tuple{Int64,UInt64}
    s2c_radix = length(boot.C2Smatrix[1])
    opermod.Qatlevel[end-s2c_radix]
end

function get_c2s_matrix(oper::CKKSOperator; scale::Real=1.0, radix::Int64=3, packlen::Int64=0)::Vector{PlainMatrix}
    if !isa(oper.packer, ComplexPackerPow2)
        throw(ErrorException("The packer must be a complex packer."))
    end

    packlen == 0 && (packlen = oper.packer.k)
    if oper.packer.k % packlen ≠ 0
        throw(DomainError("The number of slots must divide the packing parameter."))
    end

    res = Vector{PlainMatrix}(undef, radix)
    logn = trailing_zeros(packlen)
    step = round(Int64, logn / radix)
    N, ζ, group = oper.packer.N, oper.packer.ζ, oper.packer.group

    @inbounds for i = reverse(0:radix-1)
        tmpdiag = Dict{Int64,Vector{ComplexBF}}()
        tmpdiag[0] = fill(ComplexBF((scale / BigFloat(packlen, precision=192))^(1 / radix)), packlen)
        istep = 1 << (i * step)
        diagnum = min(1 << (i == radix - 1 ? logn - i * step : step), packlen ÷ istep)

        for j = 1:diagnum-1
            tmpdiag[j*istep] = zeros(ComplexBF, packlen)
            tmpdiag[packlen-j*istep] = zeros(ComplexBF, packlen)
        end

        for idx = reverse(i*step+1:(i == radix - 1 ? logn : (i + 1) * step))
            # idx-th iteration of the packing algorithm.
            len = 1 << idx
            for j1 = 0:len:packlen-1
                lenh, lenQ = len >> 1, len << 2
                gap = 2N ÷ lenQ
                for j2 = 0:lenh-1
                    idx1, idx2 = j1 + j2, j1 + j2 + lenh
                    ζi = ζ[(lenQ-(group[j2+1]&(lenQ-1)))*gap+1]

                    for idx3 = 0:istep:packlen-1
                        idx4 = (lenh + idx3) & (packlen - 1)
                        if haskey(tmpdiag, idx3) && haskey(tmpdiag, idx4)
                            t, u = tmpdiag[idx3][idx1+1], tmpdiag[idx4][idx2+1]
                            tmpdiag[idx3][idx1+1], tmpdiag[idx4][idx2+1] = t + u, (t - u) * ζi
                        end
                    end
                end
            end
        end

        # Compute the scaling factor and level.
        level = length(oper.Qatlevel) - radix + i

        evalQ = oper.operQ.evalQ
        now_Qlen, now_auxQ = oper.Qatlevel[level+1]
        tar_Qlen, tar_auxQ = oper.Qatlevel[level]

        @views nowQ = now_auxQ == 0 ? prod(evalQ.moduli[1:now_Qlen]) : prod(evalQ.moduli[1:now_Qlen-1]) * now_auxQ
        @views tarQ = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
        Δ = BigFloat(nowQ // tarQ, precision=192)

        if packlen < N / 2 && i == 0
            n1, _ = get_bsgs_param(2diagnum - 1)
            n1 *= istep
            n2 = ceil(Int64, packlen / n1)

            diag = Dict{Int64,PlainPoly}()
            for j1 = 0:n2-1
                for j2 = 0:n1-1
                    idx = j1 * n1 + j2
                    !haskey(tmpdiag, idx) && continue

                    tmpdiag[idx] = vcat(tmpdiag[idx], -im * tmpdiag[idx])
                    circshift!(tmpdiag[idx], -j1 * n1)
                    diag[idx] = encode(tmpdiag[idx], oper, level=level, scaling_factor=Δ, isntt_friendly=true)
                end
            end

            res[radix-i] = PlainMatrix(diag, (n1, n2))
        else
            n1, _ = get_bsgs_param(2diagnum - 1)
            n1 *= istep
            n2 = ceil(Int64, packlen / n1)

            diag = Dict{Int64,PlainPoly}()
            for j1 = 0:n2-1
                for j2 = 0:n1-1
                    idx = j1 * n1 + j2
                    !haskey(tmpdiag, idx) && continue

                    circshift!(tmpdiag[idx], -j1 * n1)
                    diag[idx] = encode(tmpdiag[idx], oper, level=level, scaling_factor=Δ, isntt_friendly=true)
                end
            end

            res[radix-i] = PlainMatrix(diag, (n1, n2))
        end
    end

    res
end

function get_required_key_list_c2s(radix::Int64, packlen::Int64)::Vector{Tuple{Int64}}
    logn = trailing_zeros(packlen)
    step = round(Int64, logn / radix)
    res = Vector{Tuple{Int64}}()
    @inbounds for i = reverse(0:radix-1)
        idiag = Int64[]
        istep = 1 << (i * step)
        diagnum = 1 << (i == radix - 1 ? logn - i * step : step)
        for j = 1:diagnum-1
            (j * istep) ∉ idiag && push!(idiag, j * istep)
            (packlen - j * istep) ∉ idiag && push!(idiag, packlen - j * istep)
        end

        n1, _ = get_bsgs_param(2diagnum - 1)
        n1 *= istep
        n2 = ceil(Int64, packlen / n1)

        for j = 0:n2-1
            check = false
            for k = 0:n1-1
                if (j * n1 + k) ∈ idiag
                    check = true
                    if k ≠ 0 && (k,) ∉ res
                        push!(res, (k,))
                    end
                end
            end
            if check && j ≠ 0
                push!(res, (j * n1,))
            end
        end
    end

    sort!(res, by=x -> x[1])

    res
end

function get_s2c_matrix(oper::CKKSOperator; scale::Real=1.0, radix::Int64=3, packlen::Int64=0)::Vector{PlainMatrix}
    if !isa(oper.packer, ComplexPackerPow2)
        throw(ErrorException("The packer must be a complex packer."))
    end

    packlen == 0 && (packlen = oper.packer.k)
    if oper.packer.k % packlen ≠ 0
        throw(DomainError("The number of slots must divide the packing parameter."))
    end

    res = Vector{PlainMatrix}(undef, radix)
    logn = trailing_zeros(packlen)
    step = round(Int64, logn / radix)
    N, ζ, group = oper.packer.N, oper.packer.ζ, oper.packer.group

    @inbounds for i = 0:radix-1
        tmpdiag = Dict{Int64,Vector{ComplexBF}}()
        tmpdiag[0] = fill(ComplexBF(scale)^(1 / radix), packlen)
        istep = 1 << (i * step)
        diagnum = min(1 << (i == radix - 1 ? logn - i * step : step), packlen ÷ istep)

        for j = 1:diagnum-1
            tmpdiag[j*istep] = zeros(ComplexBF, packlen)
            tmpdiag[packlen-j*istep] = zeros(ComplexBF, packlen)
        end

        for idx = i*step+1:(i == radix - 1 ? logn : (i + 1) * step)
            # idx-th iteration of the unpacking algorithm.
            len = 1 << idx
            for j1 = 0:len:packlen-1
                lenh, lenQ = len >> 1, len << 2
                gap = 2N ÷ lenQ
                for j2 = 0:lenh-1
                    idx1, idx2 = j1 + j2, j1 + j2 + lenh
                    ζi = ζ[(group[j2+1]&(lenQ-1))*gap+1]

                    for idx3 = 0:istep:packlen-1
                        idx4 = (lenh + idx3) & (packlen - 1)
                        if haskey(tmpdiag, idx3) && haskey(tmpdiag, idx4)
                            t, u = tmpdiag[idx3][idx1+1], tmpdiag[idx4][idx2+1] * ζi
                            tmpdiag[idx3][idx1+1], tmpdiag[idx4][idx2+1] = t + u, t - u
                        end
                    end
                end
            end
        end

        # Compute the scaling factor and level.
        level = length(oper.Qatlevel) - i - 1

        evalQ = oper.operQ.evalQ
        now_Qlen, now_auxQ = oper.Qatlevel[level+1]
        tar_Qlen, tar_auxQ = oper.Qatlevel[level]

        @views nowQ = now_auxQ == 0 ? prod(evalQ.moduli[1:now_Qlen]) : prod(evalQ.moduli[1:now_Qlen-1]) * now_auxQ
        @views tarQ = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
        Δ = BigFloat(nowQ // tarQ, precision=192)

        if packlen < N / 2 && i == 0
            n1, _ = get_bsgs_param(2diagnum - 1)
            n1 *= istep
            n2 = ceil(Int64, packlen / n1)

            diag = Dict{Int64,PlainPoly}()
            for j1 = 0:n2-1
                for j2 = 0:n1-1
                    idx = j1 * n1 + j2
                    !haskey(tmpdiag, idx) && continue

                    tmpdiag[idx] = vcat(tmpdiag[idx], tmpdiag[idx])

                    @. @views tmpdiag[idx][packlen+idx+1:end] *= im
                    @. @views tmpdiag[idx][1:idx] *= im

                    circshift!(tmpdiag[idx], -j1 * n1)
                    diag[idx] = encode(tmpdiag[idx], oper, level=level, scaling_factor=Δ, isntt_friendly=true)
                end
            end

            res[i+1] = PlainMatrix(diag, (n1, n2))
        else
            n1, _ = get_bsgs_param(2diagnum - 1)
            n1 *= istep
            n2 = ceil(Int64, packlen / n1)

            diag = Dict{Int64,PlainPoly}()
            for j1 = 0:n2-1
                for j2 = 0:n1-1
                    idx = j1 * n1 + j2
                    !haskey(tmpdiag, idx) && continue

                    circshift!(tmpdiag[idx], -j1 * n1)
                    diag[idx] = encode(tmpdiag[idx], oper, level=level, scaling_factor=Δ, isntt_friendly=true)
                end
            end

            res[i+1] = PlainMatrix(diag, (n1, n2))
        end
    end

    res
end

function get_required_key_list_s2c(radix::Int64, packlen::Int64)::Vector{Tuple{Int64}}
    logn = trailing_zeros(packlen)
    step = round(Int64, logn / radix)
    res = Vector{Tuple{Int64}}()
    @inbounds for i = 0:radix-1
        idiag = Int64[]
        istep = 1 << (i * step)
        diagnum = 1 << (i == radix - 1 ? logn - i * step : step)
        for j = 1:diagnum-1
            (j * istep) ∉ idiag && push!(idiag, j * istep)
            (packlen - j * istep) ∉ idiag && push!(idiag, packlen - j * istep)
        end

        n1, _ = get_bsgs_param(2diagnum - 1)
        n1 *= istep
        n2 = ceil(Int64, packlen / n1)

        for j = 0:n2-1
            check = false
            for k = 0:n1-1
                if (j * n1 + k) ∈ idiag
                    check = true
                    if k ≠ 0 && (k,) ∉ res
                        push!(res, (k,))
                    end
                end
            end
            if check && j ≠ 0
                push!(res, (j * n1,))
            end
        end
    end

    sort!(res, by=x -> x[1])

    res
end

function modraise(ct::CKKS, oper::CKKSOperator, boot::CKKSBootstrapper)::Tuple{CKKS,BigFloat}
    gap, operQ0, operc2s, ksktosparse, kskfromsparse = boot.param.gap, boot.operQ0, boot.operc2s, boot.ksktosparse, boot.kskfromsparse

    # Define the temporary ciphertext.
    evalQ = oper.operQ.evalQ
    Qlen, auxQ = length(ct), ct.val.auxQ[]
    tmp, res = oper.ct_buff[1][1:Qlen], similar(ct)

    # Scale the input ciphertext.
    @views Q0 = auxQ == 0 ? evalQ.moduli[1:Qlen] : vcat(evalQ.moduli[1:Qlen-1], Modulus(auxQ))

    # Scale the ciphertext modulus to bootQ0.
    Qdiff_real = prod(Q0) / BigFloat(gap * ct.scale[], precision=192)
    Qdiff_int = round(BigInt, Qdiff_real)

    for i = 1:Qlen
        Bmul_to!(tmp.val.b[i], (Qdiff_int % Q0[i].Q) % UInt64, ct.val.b[i], Q0[i])
        Bmul_to!(tmp.val.a[i], (Qdiff_int % Q0[i].Q) % UInt64, ct.val.a[i], Q0[i])
    end
    tmp.val.auxQ[] = ct.val.auxQ[]
    tmp.val.b.isntt[] = ct.val.b.isntt[]
    tmp.val.a.isntt[] = ct.val.a.isntt[]

    bootQ0 = operQ0.evalQ.moduli[1]
    scale_to!(tmp.val, tmp.val, 1, oper.operQ, auxQ=bootQ0.Q, isntt=false)
    tmp.scale[] = Qdiff_int / Qdiff_real * bootQ0.Q / gap

    # Sparse key encapsulation.
    tmp.val.auxQ[] = 0
    if !isnothing(ksktosparse)
        keyswitch_to!(tmp.val, tmp.val, ksktosparse, operQ0)
        intt_to!(tmp.val, tmp.val, operQ0)
    end

    # Resize the output ciphertext.
    maxlevel = length(operc2s.Qatlevel) - 1
    Qlen, auxQ = operc2s.Qatlevel[end]
    resize!(res.val, Qlen)
    res.val.auxQ[] = auxQ

    # Modulus raise.
    QL = auxQ == 0 ? oper.operQ.evalQ.moduli[1:Qlen] : vcat(oper.operQ.evalQ.moduli[1:Qlen-1], Modulus(auxQ))
    be = BasisExtender([bootQ0], QL)
    basis_extend_to!(res.val.b.coeffs, tmp.val.b.coeffs, be)
    basis_extend_to!(res.val.a.coeffs, tmp.val.a.coeffs, be)
    res.val.b.isntt[] = false
    res.val.a.isntt[] = false

    res.scale[] = bootQ0.Q
    res.level[] = maxlevel

    # Secret key encapsulation.
    if !isnothing(kskfromsparse)
        keyswitch_to!(res.val, res.val, kskfromsparse, oper.operQ)
    else
        ntt_to!(res.val, res.val, oper.operQ)
    end

    res, Qdiff_int / Qdiff_real
end

function subsum!(ct::CKKS, oper::CKKSOperator, atk::Dict{Int64,RLEV}, subringdegree::Int64)::Nothing
    N, k = oper.operQ.param.N, subringdegree

    N == k && return nothing

    logN, logk = trailing_zeros(N), trailing_zeros(k)

    Qlen, auxQ = length(ct.val), ct.val.auxQ[]
    Q = auxQ == 0 ? oper.operQ.evalQ.moduli[1:Qlen] : vcat(oper.operQ.evalQ.moduli[1:Qlen], auxQ)

    tmp = oper.ct_buff[2][1:Qlen]

    for i = 1:Qlen
        invmodQ = invmod(1 << (logN - logk - 1), Q[i])
        Bmul_to!(ct.val.b[i], invmodQ, ct.val.b[i], Q[i])
        Bmul_to!(ct.val.a[i], invmodQ, ct.val.a[i], Q[i])
    end

    for i = logN:-1:logk+2
        idx = (1 << i) + 1
        automorphism_to!(tmp, ct, idx, atk[idx], oper)
        add_to!(ct, ct, tmp, oper)
    end

    return nothing
end

function coeffs2slots(ct::CKKS, atk::Dict{Int64,RLEV}, boot::CKKSBootstrapper)::CKKS
    operc2s, C2Smatrix = boot.operc2s, boot.C2Smatrix

    res = mul(C2Smatrix[1], ct, atk, operc2s)
    for i = 2:length(C2Smatrix)
        mul_to!(res, C2Smatrix[i], res, atk, operc2s)
    end
    res.scale[] = ct.scale[] * 2^boot.param.mod_bits / boot.operQ0.evalQ.moduli[1].Q

    res
end

function evalmod(ct::CKKS, rlk::RLEV, boot::CKKSBootstrapper; scale::BigFloat=BigFloat(1.0, precision=192))::CKKS
    opermod, sinepoly, invpoly = boot.opermod, boot.sinepoly, boot.invpoly

    ct.level[] = length(opermod.Qatlevel) - 1

    sinedeg = degree(sinepoly)
    maxlevel = trailing_zeros(Base._nextpow2(sinedeg + 1))
    maxdeg = (1 << maxlevel) - 1
    Δ_target = opermod.scaling_factor

    sinebasis = compute_basis(sinepoly, ct, rlk, opermod)
    sineres = evalrecurse(0, maxdeg, sinepoly, sinebasis, Δ_target, rlk, opermod)
    if isnothing(sineres)
        throw(ErrorException("The polynomial evaluation failed at sine evaluation."))
    end

    invdeg = degree(invpoly)
    maxlevel = trailing_zeros(Base._nextpow2(invdeg + 1))
    maxdeg = (1 << maxlevel) - 1
    Δ_target = opermod.scaling_factor * scale

    invbasis = compute_basis(invpoly, sineres, rlk, opermod)
    invres = evalrecurse(0, maxdeg, invpoly, invbasis, Δ_target, rlk, opermod)
    if isnothing(invres)
        throw(ErrorException("The polynomial evaluation failed at inverse evaluation."))
    end

    return invres
end

function slots2coeffs(ct::CKKS, atk::Dict{Int64,RLEV}, boot::CKKSBootstrapper)::CKKS
    opers2c, S2Cmatrix = boot.opers2c, boot.S2Cmatrix

    ct.level[] = length(opers2c.Qatlevel) - 1
    res = mul(S2Cmatrix[1], ct, atk, opers2c)
    for i = 2:length(S2Cmatrix)
        mul_to!(res, S2Cmatrix[i], res, atk, opers2c)
    end

    res
end

function set_modulus(ct::CKKS, oper::CKKSOperator, boot::CKKSBootstrapper)::Tuple{CKKS,BigFloat}
    gap, operQ0, opers2c = boot.param.gap, boot.operQ0, boot.opers2c

    # Scale the input ciphertext.
    evalQ = oper.operQ.evalQ
    now_Qlen, now_auxQ = length(ct), ct.val.auxQ[]
    tar_Qlen, tar_auxQ = opers2c.Qatlevel[end]

    @views nowQ = now_auxQ == 0 ? evalQ.moduli[1:now_Qlen] : vcat(evalQ.moduli[1:now_Qlen-1], Modulus(now_auxQ))
    @views tarQ = tar_auxQ == 0 ? evalQ.moduli[1:tar_Qlen] : vcat(evalQ.moduli[1:tar_Qlen-1], Modulus(tar_auxQ))

    res = similar(ct)

    # Scale the ciphertext modulus to target Q.
    bootQ0 = operQ0.evalQ.moduli[1]
    Qdiff_real = prod(nowQ) / prod(tarQ) * bootQ0.Q / BigFloat(gap * ct.scale[], precision=192)
    Qdiff_int = round(BigInt, Qdiff_real)

    for i = 1:now_Qlen
        Bmul_to!(res.val.b[i], (Qdiff_int % nowQ[i].Q) % UInt64, ct.val.b[i], nowQ[i])
        Bmul_to!(res.val.a[i], (Qdiff_int % nowQ[i].Q) % UInt64, ct.val.a[i], nowQ[i])
    end
    res.val.auxQ[] = now_auxQ
    res.val.b.isntt[] = ct.val.b.isntt[]
    res.val.a.isntt[] = ct.val.a.isntt[]
    scale_to!(res.val, res.val, tar_Qlen, oper.operQ, auxQ=tar_auxQ)
    res.scale[] = Qdiff_int / Qdiff_real * bootQ0.Q / gap

    res, Qdiff_int / Qdiff_real
end

function bootstrap(ct::CKKS, oper::CKKSOperator, atk::Dict{Int64,RLEV}, rlk::RLEV, boot::CKKSBootstrapper)::CKKS
    if boot.param.isthin
        # Drop the level.
        ct_in, scale = set_modulus(ct, oper, boot)
        ct_in.level[] = length(boot.opers2c.Qatlevel) - 1

        # Slots2coeffs.
        ct_s2c = slots2coeffs(ct_in, atk, boot)

        # Modraise.
        ct_modraise, scale2 = modraise(ct_s2c, oper, boot)
        subsum!(ct_modraise, boot.operc2s, atk, boot.param.packlen)

        # Coeffs2slots.
        ct_c2s = coeffs2slots(ct_modraise, atk, boot)
        tmp = automorphism(ct_c2s, -1, atk[-1], boot.operc2s)
        add_to!(ct_c2s, ct_c2s, tmp, boot.operc2s)

        # Evalmod.
        correction = ct.scale[] / 2^boot.param.mod_bits * boot.param.gap / (scale * scale2)
        ct_evalmod = evalmod(ct_c2s, rlk, boot, scale=correction)
        ct_evalmod.scale[] = ct.scale[]

        return ct_evalmod
    else
        # Drop the level.
        ct_Q0 = drop_level(ct, 0, oper)

        # Modraise.
        ct_modraise, scale = modraise(ct_Q0, oper, boot)
        subsum!(ct_modraise, boot.operc2s, atk, boot.param.packlen)

        # Coeffs2slots.
        ct_c2s = coeffs2slots(ct_modraise, atk, boot)

        ct_c2s2 = automorphism(ct_c2s, -1, atk[-1], boot.operc2s)
        ct_c2s1 = add(ct_c2s, ct_c2s2, boot.operc2s)
        sub_to!(ct_c2s2, ct_c2s2, ct_c2s, boot.operc2s)
        mul_by_im_to!(ct_c2s2, ct_c2s2, boot.operc2s)

        # Evalmod.
        correction = ct.scale[] / 2^boot.param.mod_bits * boot.param.gap / scale
        ct_evalmod1 = evalmod(ct_c2s1, rlk, boot, scale=correction)
        ct_evalmod2 = evalmod(ct_c2s2, rlk, boot, scale=correction)

        # Slots2coeffs.
        mul_by_im_to!(ct_evalmod2, ct_evalmod2, boot.opermod)
        add_to!(ct_evalmod1, ct_evalmod1, ct_evalmod2, boot.opermod)
        ct_s2c = slots2coeffs(ct_evalmod1, atk, boot)

        # Set the scale.
        ct_s2c.scale[] = ct.scale[]

        return ct_s2c
    end
end