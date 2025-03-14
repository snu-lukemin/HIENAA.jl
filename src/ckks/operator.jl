"""
CKKSOperator is a struct for arithmetic operations over CKKS ciphertexts.
"""
struct CKKSOperator <: HEOperator
    scaling_factor::BigFloat
    operQ::Operator
    packer::ComplexPacker
    pack_buff::Vector{Int128}
    ct_buff::Vector{CKKS}
    tensor_buff::Tensor
    Qatlevel::Vector{Tuple{Int64,UInt64}}

    function CKKSOperator(param::CKKSParameters)::CKKSOperator
        ring_param, P, Q, dlen, Δ = param.ring_param, param.P, param.Q, param.dlen, BigFloat(param.scaling_factor, precision=192)

        # define operator.
        rlwe_param = RLWEParameters(ring_param, P, Q, dlen)
        operQ = Operator(rlwe_param)

        # define packer
        packer = ComplexPacker(ring_param)

        # define buff.
        pack_buff = Vector{Int128}(undef, ring_param.N)
        ct_buff = CKKS[CKKS(ring_param.N, length(Q), 0) for _ = 1:4]
        tensor_buff = Tensor(ring_param.N, length(Q) + length(P))

        # Compute the length of modulus chain, and auxQ at each level.
        minQ = log2(Δ)
        Qatlevel = Vector{Tuple{Int64,UInt64}}(undef, 0)

        Qlen = length(Q)
        auxQ = UInt64(0)
        setprecision(192)
        while Qlen > 0
            pushfirst!(Qatlevel, (Qlen, auxQ))

            nowQ = 0.0
            @inbounds for i = 1:Qlen-1
                nowQ += log2(Q[i])
            end
            nowQ += auxQ == 0 ? log2(Q[Qlen]) : log2(auxQ)
            nowQ < minQ && break

            auxQ == 0 && (auxQ = Q[Qlen])
            if Δ > auxQ
                tmpauxQ = auxQ / Δ
                while tmpauxQ < 1
                    Qlen -= 1
                    tmpauxQ *= Q[Qlen]
                end
                auxQ = round(UInt64, tmpauxQ)
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            else
                auxQ = round(UInt64, 1 / Δ * auxQ)
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end
        popfirst!(Qatlevel)

        new(Δ, operQ, packer, pack_buff, ct_buff, tensor_buff, Qatlevel)
    end

    # Used for plaintext switching.
    function CKKSOperator(oper::CKKSOperator, scaling_factor::Real, top_Qlen::Int64, top_auxQ::UInt64)::CKKSOperator
        # Define others.
        operQ, pack_buff, ct_buff, tensor_buff = oper.operQ, oper.pack_buff, oper.ct_buff, oper.tensor_buff
        packer = oper.packer
        Δ = BigFloat(scaling_factor, precision=192)

        # Compute the length of modulus chain, and auxQ at each level.
        minQ = log2(Δ)
        Qatlevel = Vector{Tuple{Int64,UInt64}}(undef, 0)

        Q = operQ.evalQ.moduli
        Qlen, auxQ = top_Qlen, top_auxQ
        while true
            pushfirst!(Qatlevel, (Qlen, auxQ))

            nowQ = 0.0
            @inbounds for i = 1:Qlen-1
                nowQ += log2(Q[i].Q)
            end
            nowQ += auxQ == 0 ? log2(Q[Qlen].Q) : log2(auxQ)
            nowQ < minQ && break

            auxQ == 0 && (auxQ = Q[Qlen].Q)
            if Δ > auxQ
                tmpauxQ = auxQ / Δ
                while tmpauxQ < 1
                    Qlen -= 1
                    tmpauxQ *= Q[Qlen].Q
                end
                auxQ = round(UInt64, tmpauxQ)
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            else
                auxQ = round(UInt64, 1 / Δ * auxQ)
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end
        popfirst!(Qatlevel)

        new(Δ, operQ, packer, pack_buff, ct_buff, tensor_buff, Qatlevel)
    end
end

#======================================================================================#
################################## ENCODE AND DECODE ###################################
#======================================================================================#

function encode(msg::Real, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false, scaling_factor::Real=0.0)::PlainConst
    if level < 0
        throw(DomainError("The level should be greater than 0."))
    end

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen, _ = Qatlevel[level+1]

    if isPQ
        if ismissing(oper.operQ.evalP)
            throw(ErrorException("The parameter does not support PQ."))
        end
        res = PlainConst(Qlen + length(oper.operQ.evalP))
        encode_to!(res, msg, oper, level=level, isPQ=true, scaling_factor=scaling_factor)
    else
        res = PlainConst(Qlen)
        encode_to!(res, msg, oper, level=level, isPQ=false, scaling_factor=scaling_factor)
    end

    res
end

function encode_to!(res::PlainConst, msg::Real, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false, scaling_factor::Real=0.0)::Nothing
    if level < 0
        throw(DomainError("The level should be greater than 0."))
    end

    if scaling_factor == 0.0
        Δ = oper.scaling_factor
    else
        Δ = scaling_factor
    end

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tmp = round(Int128, msg * Δ)

    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    if isPQ
        if ismissing(operQ.evalP)
            throw(ErrorException("The parameter does not support PQ."))
        end

        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ, isPQ=true)
        Plen = length(operQ.evalP)

        resize!(res, Plen + Qlen)
        for i = 1:Plen+Qlen
            res.val.vals[i] = Bred(tmp, evalQ[i])
        end
        res.auxQ[] = auxQ
        res.isPQ[] = true
        res.scale[] = Δ
    else
        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)

        resize!(res, Qlen)
        for i = 1:Qlen
            res.val.vals[i] = Bred(tmp, evalQ[i])
        end
        res.auxQ[] = auxQ
        res.isPQ[] = false
        res.scale[] = Δ
    end

    return nothing
end

function encode(msg::Vector{<:Union{Real,Complex{<:Real}}}, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true, scaling_factor::Real=0.0, isntt_friendly::Bool=false)::PlainPoly
    if level < 0
        throw(DomainError("The level should be greater than 0."))
    end

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen, _ = Qatlevel[level+1]

    if isPQ
        if ismissing(oper.operQ.evalP)
            throw(ErrorException("The parameter does not support PQ."))
        end
        res = PlainPoly(oper.operQ.param.N, Qlen + length(oper.operQ.evalP))
        encode_to!(res, msg, oper, level=level, isPQ=true, ispacking=ispacking, scaling_factor=scaling_factor, isntt_friendly=isntt_friendly)
    else
        res = PlainPoly(oper.operQ.param.N, Qlen)
        encode_to!(res, msg, oper, level=level, isPQ=false, ispacking=ispacking, scaling_factor=scaling_factor, isntt_friendly=isntt_friendly)
    end

    res
end

function encode_to!(res::PlainPoly, msg::Vector{<:Union{Real,Complex{<:Real}}}, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true, scaling_factor::Real=0.0, isntt_friendly::Bool=false)::Nothing
    if level < 0
        throw(DomainError("The level should be greater than 0."))
    end

    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    if scaling_factor == 0.0
        Δ = oper.scaling_factor
    else
        Δ = scaling_factor
    end

    if ispacking
        packer, pack_buff = oper.packer, oper.pack_buff
        pack_to!(pack_buff, msg, Δ, packer)
    else
        if length(pack_buff) ≠ length(msg)
            throw(DimensionMismatch("The length of the message should match the ring size."))
        end
        @. pack_buff = round(Int128, msg * Δ)
    end

    if isPQ
        if ismissing(operQ.evalP)
            throw(ErrorException("The parameter does not support PQ."))
        end

        Qlen, auxQ = Qatlevel[level+1]
        isntt_friendly && (auxQ = UInt64(0))
        Plen = length(operQ.evalP)
        evalPQ = geteval_at(Plen + Qlen, operQ, auxQ=auxQ, isPQ=true)

        resize!(res, Plen + Qlen)
        for i = eachindex(evalPQ)
            Bred_to!(res.val.coeffs[i], pack_buff, evalPQ[i].Q)
        end

        res.val.isntt[] = false
        ntt!(res.val, evalPQ)

        res.auxQ[] = auxQ
        res.isPQ[] = true
        res.scale[] = Δ
    else
        Qlen, auxQ = Qatlevel[level+1]
        isntt_friendly && (auxQ = UInt64(0))
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)

        resize!(res, Qlen)
        for i = eachindex(evalQ)
            Bred_to!(res.val.coeffs[i], pack_buff, evalQ[i].Q)
        end

        res.val.isntt[] = false
        ntt!(res.val, evalQ)

        res.auxQ[] = auxQ
        res.isPQ[] = false
        res.scale[] = Δ
    end

    return nothing
end

#TODO decode to Int128. Can it be done without BigInteger arithmetic?
function decode(x::PlainConst, oper::CKKSOperator)::BigFloat
    to_big(x, oper.operQ) / x.scale[]
end

function decode(x::PlainPoly, oper::CKKSOperator; ispacking::Bool=true)::Vector{<:Union{BigFloat,ComplexBF}}
    if ispacking
        res = similar(oper.packer.buff)
    else
        res = Vector{BigFloat}(undef, x.val.N)
    end

    decode_to!(res, x, oper, ispacking=ispacking)

    res
end

function decode_to!(res::Vector{<:Union{BigFloat,ComplexBF}}, x::PlainPoly, oper::CKKSOperator; ispacking::Bool=true)::Nothing
    isPQ, auxQ = x.isPQ[], x.auxQ[]
    evalQ = geteval_at(length(x), oper.operQ, auxQ=auxQ, isPQ=isPQ)

    buff = oper.tensor_buff.vals[1][1:length(x)]
    copy!(buff, x.val)
    buff.isntt[] && intt!(buff, evalQ)

    packer, pack_buff, Δ = oper.packer, oper.pack_buff, x.scale[]
    pack_buff .= to_big(x, oper.operQ)

    if ispacking
        if length(res) ≠ oper.packer.k
            throw(DimensionMismatch("The length of the plaintext should match the ring size."))
        end
        unpack_to!(res, pack_buff, Δ, packer)
    else
        if length(res) ≠ x.val.N
            throw(DimensionMismatch("The length of the plaintext should match the ring size."))
        end
        @. res = pack_buff / Δ
    end

    return nothing
end

#======================================================================================#
################################## ENCRYPT AND DECRYPT #################################
#======================================================================================#

function encrypt(msg::PlainText, entor::Encryptor, oper::CKKSOperator)::CKKS
    res = CKKS(RLWE(oper.operQ.param.N, length(msg)), 0, 0.0)
    encrypt_to!(res, msg, entor, oper)
    res
end

function encrypt_to!(res::CKKS, msg::PlainConst, entor::Encryptor, oper::CKKSOperator)::Nothing
    Qlen, auxQ = length(msg), msg.auxQ[]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = auxQ
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # Compute (b + m, a), which has the phase m + e mod Q.
    @inbounds for i = 1:Qlen
        add_to!(val.b.coeffs[i], val.b.coeffs[i], msg.val[i], val.b.isntt[], evalQ[i])
    end

    # Set scale
    res.scale[] = msg.scale[]

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i][1] == Qlen && Qatlevel[i][2] == auxQ
            res.level[] = i - 1
            return nothing
        end
    end
    throw(ErrorException("There is no such level."))
end

function encrypt_to!(res::CKKS, msg::PlainPoly, entor::Encryptor, oper::CKKSOperator)::Nothing
    Qlen, auxQ = length(msg), msg.auxQ[]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = auxQ
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # NTT Buffer.
    buff_ntt = oper.tensor_buff[1].coeffs[1]

    # Compute (b + m, a), which has the phase m + e mod Q.
    @inbounds for i = 1:Qlen
        if !msg.val.isntt[]
            ntt_to!(buff_ntt, msg.val.coeffs[i], evalQ[i])
            add_to!(res.val.b.coeffs[i], res.val.b.coeffs[i], buff_ntt, evalQ[i])
        else
            add_to!(res.val.b.coeffs[i], res.val.b.coeffs[i], msg.val.coeffs[i], evalQ[i])
        end
    end

    # Set scale
    res.scale[] = msg.scale[]

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i][1] == Qlen && Qatlevel[i][2] == auxQ
            res.level[] = i - 1
            return nothing
        end
    end
    throw(ErrorException("There is no such level."))
end

function decrypt(ct::CKKS, entor::Encryptor, oper::CKKSOperator)::PlainPoly
    res = PlainPoly(oper.operQ.param.N, length(ct.val))
    decrypt_to!(res, ct, entor, oper)
    res
end

function decrypt_to!(res::PlainPoly, ct::CKKS, entor::Encryptor, oper::CKKSOperator)::Nothing
    resize!(res, length(ct.val))
    phase_to!(res, ct.val, entor)
    res.scale[] = ct.scale[]
    return nothing
end

#=================================================================================================#
#################################### PLAINTEXT OPERATIONS #########################################
#=================================================================================================#

#TODO plaintext addtion/subtraction/multiplication

change_level(x::PlainText, targetlvl::Integer, oper::CKKSOperator)::PlainText = begin
    res = similar(x)
    change_level_to!(res, x, targetlvl, oper)
    res
end

function change_level_to!(res::PlainPoly, x::PlainPoly, targetlvl::Integer, oper::CKKSOperator)::Nothing
    if targetlvl < 0
        throw(DomainError("The target level should be greater than 0."))
    end

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]
    change_modulus_to!(res, x, tar_Qlen, operQ, auxQ=tar_auxQ)

    return nothing
end

function change_level_to!(res::PlainConst, x::PlainConst, targetlvl::Integer, oper::CKKSOperator)::Nothing
    if targetlvl < 0
        throw(DomainError("The target level should be greater than 0."))
    end

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]
    change_modulus_to!(res, x, tar_Qlen, operQ, auxQ=tar_auxQ)

    return nothing
end

#================================================================================================#
#################################### CIPHERTEXT OPERATIONS ########################################
#================================================================================================#

drop_level(x::CKKS, targetlvl::Integer, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    drop_level_to!(res, x, targetlvl, oper)
    res
end

function drop_level_to!(res::CKKS, x::CKKS, targetlvl::Integer, oper::CKKSOperator)::Nothing
    currentlvl = x.level[]

    if targetlvl < 0
        throw(DomainError("The target level should be greater than 0."))
    end
    if targetlvl > currentlvl
        throw(DomainError("The target level should be less than or equal to the current level."))
    end

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    if now_Qlen ≠ length(x.val) || now_auxQ ≠ x.val.auxQ[]
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Easy case.
    if targetlvl == currentlvl
        resize!(res.val, length(x.val))
        copy!(res, x)
        return nothing
    end

    # Copy the ciphertext to the buffer.
    buff = operQ.ct_buff[1][1:now_Qlen]
    copy!(buff, x.val[1:now_Qlen])

    # Multiply by the divided bits.
    evalQ = geteval_at(now_Qlen, operQ, auxQ=now_auxQ)
    @views tarQ = tar_auxQ == 0 ? prod(operQ.evalQ.moduli[1:tar_Qlen]) : prod(operQ.evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
    Qdiff_real = prod(evalQ.moduli) // tarQ
    Qdiff_int = round(BigInt, Qdiff_real)
    for i = eachindex(evalQ)
        diffi = (Qdiff_int % evalQ[i].Q.Q) % UInt64
        mul_to!(buff.b[i], diffi, buff.b[i], evalQ[i])
        mul_to!(buff.a[i], diffi, buff.a[i], evalQ[i])
    end

    # Rational rescale to the new modulus.
    resize!(res.val, tar_Qlen)
    scale_to!(res.val, buff, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Update the level of the ciphertext.
    res.level[] = targetlvl

    # Update the scale of the ciphertext.
    res.scale[] = x.scale[] * Qdiff_int / Qdiff_real

    return nothing
end

rescale(x::CKKS, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    rescale_to!(res, x, oper)
    res
end

function rescale_to!(res::CKKS, x::CKKS, oper::CKKSOperator)::Nothing
    currentlvl = x.level[]

    if currentlvl <= 0
        throw(DomainError("The level of ciphertext should be greater than 0."))
    end

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[currentlvl]

    # Sanity check
    if now_Qlen ≠ length(x.val) || now_auxQ ≠ x.val.auxQ[]
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Copy the ciphertext to the buffer.
    buff = oper.ct_buff[1][1:now_Qlen]
    copy!(buff, x)

    # Rational rescale to the new modulus.
    scale_to!(res.val, buff.val, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Update the level of the ciphertext.
    res.level[] = currentlvl - 1

    # Update the scale of the ciphertext.
    evalQ = oper.operQ.evalQ
    @views nowQ = now_auxQ == 0 ? prod(evalQ.moduli[1:now_Qlen]) : prod(evalQ.moduli[1:now_Qlen-1]) * now_auxQ
    @views tarQ = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
    res.scale[] = x.scale[] * (tarQ // nowQ)

    return nothing
end

function neg(x::CKKS, oper::CKKSOperator)::CKKS
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::CKKS, x::CKKS, oper::CKKSOperator)::Nothing
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Compute res = -x.
    resize!(res.val, tar_Qlen)
    for i = 1:tar_Qlen
        neg_to!(res.val.b[i], x.val.b[i], evalQ[i])
        neg_to!(res.val.a[i], x.val.a[i], evalQ[i])
    end
    res.level[] = level
    res.scale[] = x.scale[]

    return nothing
end

add(x::CKKS, y::PlainPoly, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainPoly, y::CKKS, oper::CKKSOperator)::CKKS = add(y, x, oper)

function add_to!(res::CKKS, x::CKKS, y::PlainPoly, oper::CKKSOperator)::Nothing
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    if length(x) ≠ Qlen || x.val.auxQ[] ≠ auxQ
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Match the level.
    buff = PlainPoly(oper.tensor_buff.vals[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute x + y.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
    res.scale[] = (x.scale[] + y.scale[]) / 2

    return nothing
end

add_to!(res::CKKS, x::PlainPoly, y::CKKS, oper::CKKSOperator)::Nothing = begin
    add_to!(res, y, x, oper)
    return nothing
end

add(x::CKKS, y::PlainConst, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainConst, y::CKKS, oper::CKKSOperator)::CKKS = add(y, x, oper)

function add_to!(res::CKKS, x::CKKS, y::PlainConst, oper::CKKSOperator)::Nothing
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    if length(x) ≠ Qlen || x.val.auxQ[] ≠ auxQ
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Match the level.
    buff = change_level(y, level, oper)

    # Compute x + y.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
    res.scale[] = (x.scale[] + y.scale[]) / 2

    return nothing
end

add_to!(res::CKKS, x::PlainConst, y::CKKS, oper::CKKSOperator)::Nothing = begin
    add_to!(res, y, x, oper)
    return nothing
end

add(x::CKKS, y::CKKS, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

"""
Add two CKKS ciphertexts.
"""
function add_to!(res::CKKS, x::CKKS, y::CKKS, oper::CKKSOperator)::Nothing
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]
    tmpx, tmpy = oper.ct_buff[1][1:tar_Qlen], oper.ct_buff[2][1:tar_Qlen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    resize!(res.val, tar_Qlen)
    res.val.auxQ[] = tar_auxQ
    add_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
    res.scale[] = (x.scale[] + y.scale[]) / 2

    return nothing
end

sub(x::CKKS, y::PlainPoly, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::CKKS, x::CKKS, y::PlainPoly, oper::CKKSOperator)::Nothing
    buff = PlainPoly(oper.tensor_buff.vals[2][1:length(x)])
    neg_to!(buff, y, oper.operQ)
    add_to!(res, x, buff, oper)
    return nothing
end

sub(x::PlainPoly, y::CKKS, oper::CKKSOperator)::CKKS = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::CKKS, x::PlainPoly, y::CKKS, oper::CKKSOperator)::Nothing
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
    return nothing
end

sub(x::CKKS, y::PlainConst, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::CKKS, x::CKKS, y::PlainConst, oper::CKKSOperator)::Nothing = begin
    tmp = neg(y, oper.operQ)
    add_to!(res, x, tmp, oper)
    return nothing
end

sub(x::PlainConst, y::CKKS, oper::CKKSOperator)::CKKS = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::CKKS, x::PlainConst, y::CKKS, oper::CKKSOperator)::Nothing = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
    return nothing
end

sub(x::CKKS, y::CKKS, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

"""
Subtract two CKKS ciphertexts.
"""
function sub_to!(res::CKKS, x::CKKS, y::CKKS, oper::CKKSOperator)::Nothing
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]
    tmpx, tmpy = oper.ct_buff[1][1:tar_Qlen], oper.ct_buff[2][1:tar_Qlen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    resize!(res.val, tar_Qlen)
    res.val.auxQ[] = tar_auxQ
    sub_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
    res.scale[] = (x.scale[] + y.scale[]) / 2

    return nothing
end

mul(x::CKKS, y::PlainPoly, oper::CKKSOperator; islazy::Bool=false)::CKKS = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainPoly, y::CKKS, oper::CKKSOperator; islazy::Bool=false)::CKKS = mul(y, x, oper, islazy=islazy)

function mul_to!(res::CKKS, x::CKKS, y::PlainPoly, oper::CKKSOperator; islazy::Bool=false)::Nothing
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    if length(x) ≠ Qlen || x.val.auxQ[] ≠ auxQ
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Match the level.
    buff = PlainPoly(oper.tensor_buff.vals[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute x * y.
    resize!(res.val, Qlen)
    mul_to!(res.val, x.val, buff, operQ)
    res.level[] = level
    res.scale[] = x.scale[] * y.scale[]

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)

    return nothing
end

mul_to!(res::CKKS, x::PlainPoly, y::CKKS, oper::CKKSOperator; islazy::Bool=false)::Nothing = begin
    mul_to!(res, y, x, oper, islazy=islazy)
    return nothing
end

mul(x::CKKS, y::PlainConst, oper::CKKSOperator; islazy::Bool=false)::CKKS = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainConst, y::CKKS, oper::CKKSOperator; islazy::Bool=false)::CKKS = mul(y, x, oper, islazy=islazy)

function mul_to!(res::CKKS, x::CKKS, y::PlainConst, oper::CKKSOperator; islazy::Bool=false)::Nothing
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    if length(x) ≠ Qlen || x.val.auxQ[] ≠ auxQ
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Match the level.
    buff = change_level(y, level, oper)

    # Compute x * y.
    resize!(res.val, Qlen)
    mul_to!(res.val, x.val, buff, operQ)
    res.level[] = level
    res.scale[] = x.scale[] * y.scale[]

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)

    return nothing
end

mul_to!(res::CKKS, x::PlainConst, y::CKKS, oper::CKKSOperator; islazy::Bool=false)::Nothing = begin
    mul_to!(res, y, x, oper, islazy=islazy)
    return nothing
end

function mul(x::CKKS, y::CKKS, rlk::RLEV, oper::CKKSOperator; islazy::Bool=false)::CKKS
    res = similar(x)
    mul_to!(res, x, y, rlk, oper, islazy=islazy)
    res
end

function mul_to!(res::CKKS, x::CKKS, y::CKKS, rlk::RLEV, oper::CKKSOperator; islazy::Bool=false)::Nothing
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    if targetlvl <= 0
        throw(DomainError("The input ciphertexts should be at least at level 1."))
    end

    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Drop the unnecessary levels of input ciphertexts.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.ct_buff[1][1:xlen], oper.ct_buff[2][1:ylen]
    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    # Tensor the input ciphertexts.
    buff = oper.tensor_buff[1:3, 1:tar_Qlen]
    _tensor_to!(buff, tmpx, tmpy, oper)

    # Relinearisation.
    resize!(res.val, tar_Qlen)
    res.val.auxQ[] = tar_auxQ
    _relinearise_to!(res.val, buff, rlk, oper)
    res.level[] = targetlvl

    # Compute the final scale.
    res.scale[] = tmpx.scale[] * tmpy.scale[]

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)

    return nothing
end

function _tensor_to!(res::Tensor, x::CKKS, y::CKKS, oper::CKKSOperator)::Nothing
    xlevel, ylevel = x.level[], y.level[]
    if xlevel ≠ ylevel
        throw(DomainError("The input ciphertexts should be at the same level."))
    end

    targetlvl = xlevel
    if targetlvl <= 0
        throw(DomainError("The input ciphertexts should be at least at level 1."))
    end

    # Sanity check
    _, len = size(res)
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]
    if len ≠ tar_Qlen
        throw(DomainError("The length of the tensor should match the length of the ciphertexts."))
    end

    # define operator.
    operQ = oper.operQ

    # Tensor.
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    !x.val.b.isntt[] && ntt!(x.val.b, evalQ)
    !x.val.a.isntt[] && ntt!(x.val.a, evalQ)
    !y.val.b.isntt[] && ntt!(y.val.b, evalQ)
    !y.val.a.isntt[] && ntt!(y.val.a, evalQ)

    mul_to!(res.vals[1], x.val.b, y.val.b, evalQ)
    mul_to!(res.vals[2], x.val.a, y.val.b, evalQ)
    muladd_to!(res.vals[2], x.val.b, y.val.a, evalQ)
    mul_to!(res.vals[3], x.val.a, y.val.a, evalQ)

    res.auxQ[] = tar_auxQ

    return nothing
end

function _relinearise_to!(res::RLWE, x::Tensor, rlk::RLEV, oper::CKKSOperator)::Nothing
    # define operator.
    operQ = oper.operQ

    # Get the length of the modulus chain and the auxiliary modulus.
    _, Qlen = size(x)

    # Relinearise.
    resize!(res, Qlen)
    relinearise_to!(res, x, rlk, operQ)

    return nothing
end

function decompose_a(x::CKKS, oper::CKKSOperator)::Tensor
    len, decer = length(x.val), oper.operQ.decer

    dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.operQ.evalP)
    deca = Tensor(oper.operQ.param.N, len + Plen, dlen)

    decompose_a_to!(deca, x, oper)

    deca
end

function decompose_a_to!(deca::Tensor, x::CKKS, oper::CKKSOperator)::Nothing
    targetlvl = x.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    if length(x.val) ≠ tar_Qlen || x.val.auxQ[] ≠ tar_auxQ
        throw(DomainError("The length of the tensor should match the length of the ciphertexts."))
    end

    # Decompose.
    poly = PlainPoly(x.val.a, auxQ=tar_auxQ)
    decompose_to!(deca, poly, oper.operQ)

    return nothing
end

keyswitch(x::CKKS, ksk::RLEV, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    keyswitch_to!(res, x, ksk, oper)
    res
end

function keyswitch_to!(res::CKKS, ct::CKKS, ksk::RLEV, oper::CKKSOperator)::Nothing
    targetlvl = ct.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    if length(ct.val) ≠ tar_Qlen || ct.val.auxQ[] ≠ tar_auxQ
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Key switching.
    resize!(res.val, tar_Qlen)
    keyswitch_to!(res.val, ct.val, ksk, oper.operQ)

    res.level[] = ct.level[]
    res.scale[] = ct.scale[]

    return nothing
end

rotate(x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, oper)
    res
end::CKKS

function rotate_to!(res::CKKS, x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator)::Nothing where {N}
    packer = oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, oper)

    return nothing
end

hoisted_rotate(adec::Tensor, x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N} = begin
    res = similar(x)
    hoisted_rotate_to!(res, adec, x, idx, rtk, oper)
    res
end::CKKS

hoisted_rotate_to!(res::CKKS, adec::Tensor, x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N} = begin
    packer = oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, x, autidx, rtk, oper)

    return nothing
end::Nothing

#TODO conjugation.

automorphism(x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

function automorphism_to!(res::CKKS, x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator)::Nothing
    targetlvl = x.level[]

    # Get the length of the modulus chain.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    if length(x.val) ≠ tar_Qlen || x.val.auxQ[] ≠ tar_auxQ
        throw(DomainError("The length of the tensor should match the length of the ciphertexts."))
    end

    # Automorphism.
    resize!(res.val, tar_Qlen)
    automorphism_to!(res.val, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
    res.scale[] = x.scale[]

    return nothing
end

hoisted_automorphism(adec::Tensor, x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    hoisted_automorphism_to!(res, adec, x, idx, atk, oper)
    res
end

function hoisted_automorphism_to!(res::CKKS, adec::Tensor, x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator)::Nothing
    targetlvl = x.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    if length(x.val) ≠ tar_Qlen || x.val.auxQ[] ≠ tar_auxQ
        throw(DomainError("The length of the tensor should match the length of the ciphertexts."))
    end

    # Automorphism.
    resize!(res.val, tar_Qlen)
    hoisted_automorphism_to!(res.val, adec, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
    res.scale[] = x.scale[]

    return nothing
end

# Used for bootstrapping.
function mul_by_im_to!(res::CKKS, x::CKKS, oper::CKKSOperator)::Nothing
    if !ispow2(oper.operQ.param.N)
        throw(DomainError("The ring size should be even."))
    end

    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Compute res = x * i.
    N = oper.packer.N
    halfN = N >> 1
    resize!(res.val, tar_Qlen)
    copy!(res, x)
    for i = 1:tar_Qlen
        evalQi = evalQ[i]
        if res.val.b.isntt[] && res.val.a.isntt[] && typeof(evalQi) == PolyEvaluatorNTT
            ζ = iMform(evalQi.ntter.Ψ[halfN+1], evalQi.Q)
            imoverQ = powermod(ζ, halfN, evalQi.Q)

            @views mul_to!(res.val.b[i][1:halfN], imoverQ, res.val.b[i][1:halfN], evalQi)
            @views mul_to!(res.val.b[i][halfN+1:end], neg(imoverQ, evalQi), res.val.b[i][halfN+1:end], evalQi)

            @views mul_to!(res.val.a[i][1:halfN], imoverQ, res.val.a[i][1:halfN], evalQi)
            @views mul_to!(res.val.a[i][halfN+1:end], neg(imoverQ, evalQi), res.val.a[i][halfN+1:end], evalQi)
        else
            circshift!(res.val.b[i], halfN)
            @views neg_to!(res.val.b[i][1:halfN], res.val.b[i][1:halfN], evalQi)

            circshift!(res.val.a[i], halfN)
            @views neg_to!(res.val.a[i][1:halfN], res.val.a[i][1:halfN], evalQi)
        end
    end

    res.level[] = level
    res.scale[] = x.scale[]

    return nothing
end

#================================================================================#
##################################### Matrix #####################################
#================================================================================#

(::Type{PlainMatrix})(M::Matrix{<:Complex{<:AbstractFloat}}, level::Integer, oper::CKKSOperator; BSGSparam::Tuple{Int64,Int64}=(0, 0))::PlainMatrix = begin
    row, col = size(M)
    if row ≠ col
        throw(DomainError("Currently, only square matrices are supported."))
    end
    if ismissing(oper.packer)
        throw(ErrorException("The operator must have a packer."))
    end
    if !isa(oper.packer, ComplexPackerPow2)
        throw(ErrorException("The packer must be a subring packer."))
    end
    if oper.packer.k % row ≠ 0
        throw(DomainError("The number of rows must divide the packing parameter."))
    end

    res = Dict{Int64,PlainPoly}()
    mattmp = Vector{ComplexBF}(undef, row)
    n1, n2 = BSGSparam == (0, 0) ? get_bsgs_param(row) : BSGSparam
    if n1 * n2 < row
        throw(DomainError("The BSGS parameters are too small."))
    end

    # Set the scaling factor.
    evalQ = oper.operQ.evalQ
    now_Qlen, now_auxQ = oper.Qatlevel[level+1]
    tar_Qlen, tar_auxQ = oper.Qatlevel[level]

    @views nowQ = now_auxQ == 0 ? prod(evalQ.moduli[1:now_Qlen]) : prod(evalQ.moduli[1:now_Qlen-1]) * now_auxQ
    @views tarQ = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
    Δ = BigFloat(nowQ // tarQ, precision=192)

    # Encode.
    packer, pack_buff = oper.packer, oper.pack_buff
    for i = 0:n2-1
        for j = 0:n1-1
            i * n1 + j ≥ row && break

            @inbounds for k = 0:row-1
                mattmp[k+1] = M[k+1, mod(k - i * n1 - j, row)+1]
            end
            all(x -> isapprox(x, 0, atol=1 / Δ), mattmp) && continue

            circshift!(mattmp, -i * n1)
            res[i*n1+j] = encode(mattmp, oper, level=level, scaling_factor=Δ, isntt_friendly=true)
        end
    end

    PlainMatrix(res, (n1, n2))
end

mul(M::PlainMatrix, x::CKKS, atk::Dict{Int64,RLEV}, oper::CKKSOperator)::CKKS = begin
    res = similar(x)
    mul_to!(res, M, x, atk, oper)
    res
end

mul_to!(res::CKKS, M::PlainMatrix, x::CKKS, atk::Dict{Int64,RLEV}, oper::CKKSOperator)::Nothing = begin
    # Sanity check.
    xlevel = x.level[]
    if xlevel <= 0
        throw(DomainError("The input ciphertext should be at least at level 1."))
    end

    now_Qlen, now_auxQ = oper.Qatlevel[xlevel+1]
    tar_Qlen, tar_auxQ = oper.Qatlevel[xlevel]
    if now_Qlen ≠ length(x) || now_auxQ ≠ x.val.auxQ[]
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Copy the value.
    resize!(res.val, now_Qlen)
    copy!(res, x)

    # Scale the ciphertext to NTT-friendly modulus.
    operQ = oper.operQ
    scale_to!(res.val, res.val, now_Qlen, operQ)

    # Perform the matrix multiplication.
    _mul_RLWE!(M, res.val, atk, oper)

    # Rescale the ciphertext.
    scale_to!(res.val, res.val, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Update the level.
    res.level[] = xlevel - 1

    return nothing
end

#================================================================================#
############################# Polynomial Evaluation ##############################
#================================================================================#

(::Type{PolyCoeffs})(coeffs::Vector{<:AbstractFloat}, oper::CKKSOperator; type::Symbol=:monomial)::PolyCoeffs = begin
    if type ≠ :monomial && type ≠ :chebyshev
        throw(DomainError("Unknown polynomial types."))
    end

    setprecision(192)
    if type == :monomial
        res = Dict{Int64,BigFloat}()
        for i = eachindex(coeffs)
            if abs(coeffs[i]) ≥ 1 / oper.scaling_factor
                res[i-1] = coeffs[i]
            end
        end
    elseif type == :chebyshev
        deg = length(coeffs) - 1
        maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
        babydeg = 1 << ceil(Int64, maxlevel / 2)
        babylevel = trailing_zeros(babydeg)
        tmpcoeffs = vcat(coeffs, zeros(typeof(coeffs[1]), 1 << maxlevel - deg - 1))

        @views for i = reverse(1:maxlevel-1)
            d = 1 << i

            if i ≥ babylevel
                for j = 0:1<<(maxlevel-i-1)-1
                    hi = 2d * (j + 1)
                    @. tmpcoeffs[hi-d:-1:hi-2d+2] -= tmpcoeffs[hi-d+2:hi]
                    @. tmpcoeffs[hi-d+2:hi] *= 2
                end
            else
                @. tmpcoeffs[end-d:-1:end-2d+2] -= tmpcoeffs[end-d+2:end]
                @. tmpcoeffs[end-d+2:end] *= 2
            end
        end
        resize!(tmpcoeffs, deg + 1)

        res = Dict{Int64,BigFloat}()
        for i = eachindex(tmpcoeffs)
            if abs(tmpcoeffs[i]) ≥ 1 / oper.scaling_factor
                res[i-1] = tmpcoeffs[i]
            end
        end
    end

    PolyCoeffs(res, type)
end

function evaluate(poly::PolyCoeffs, x::CKKS, rlk::RLEV, oper::CKKSOperator)::CKKS
    if poly.type ≠ :monomial && poly.type ≠ :chebyshev
        throw(DomainError("Unknown polynomial types."))
    end

    # Compute the basis for babystep.
    basis = compute_basis(poly, x, rlk, oper)

    # Evaluate the polynomial using BSGS algorithm.
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    maxdeg = (1 << maxlevel) - 1
    Δ_target = BigFloat(oper.scaling_factor, precision=192)
    res = evalrecurse(0, maxdeg, poly, basis, Δ_target, rlk, oper)
    
    if isnothing(res)
        throw(ErrorException("The polynomial evaluation failed."))
    end

    return res
end

function compute_basis(poly::PolyCoeffs, x::CKKS, rlk::RLEV, oper::CKKSOperator)::Dict{Int64,CKKS}
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    maxdeg = (1 << maxlevel) - 1
    babydeg = 1 << ceil(Int64, maxlevel / 2)
    babylevel = trailing_zeros(babydeg)
    basis = Dict{Int64,CKKS}()

    if poly.type == :monomial
        # Compute the power-of-two monomials.
        basis[1] = deepcopy(x)
        for i = 1:maxlevel-1
            basis[1<<i] = mul(basis[1<<(i-1)], basis[1<<(i-1)], rlk, oper)
        end

        @inbounds for i = 1:babylevel-1
            # Compute the power-of-two monomial.
            for j = 1:(1<<i)-1
                # Checks if the polynomial is sparse.
                checkj = !isassigned(poly, (1 << i) + j)

                # Used in baby step generation?
                for idx = i+1:babylevel-1
                    checkj = checkj && !isassigned(poly, (1 << idx) + (1 << i) + j)
                end
                # Used in the final evaluation?
                for idx = 0:maxdeg÷babydeg-2
                    checkj = checkj && !isassigned(poly, idx * babydeg + j)
                end
                checkj && continue

                # Compute the monomial if needed in the future.
                basis[(1<<i)+j] = mul(basis[j], basis[1<<i], rlk, oper)
            end
        end
    elseif poly.type == :chebyshev
        # Precompute constants.
        xlen = length(x)
        two, one = PlainConst(xlen), PlainConst(xlen)

        # Compute the power-of-two bases.
        # T_2^i(x) = 2 * T_2^{i-1}(x) * T_2^{i-1}(x) - 1
        basis[1] = x
        for i = 1:maxlevel-1
            basis[1<<i] = mul(basis[1<<(i-1)], basis[1<<(i-1)], rlk, oper)
            encode_to!(two, 2, oper, level=basis[1<<i].level[], scaling_factor=1.0)
            encode_to!(one, 1, oper, level=basis[1<<i].level[], scaling_factor=basis[1<<i].scale[])
            mul_to!(basis[1<<i], basis[1<<i], two, oper, islazy=true)
            sub_to!(basis[1<<i], basis[1<<i], one, oper)
        end

        operQ, evalQ = oper.operQ, oper.operQ.evalQ
        @inbounds for i = 1:babylevel-1
            for j = 1:(1<<i)-1
                # Checks if the polynomial is sparse.
                checkj = !isassigned(poly, (1 << i) + j)

                # Used in baby step generation?
                for idx = i+1:babylevel-1
                    checkj = checkj && !isassigned(poly, (1 << idx) + (1 << i) + j)
                end
                # Used in the final evaluation?
                for idx = 0:maxdeg÷babydeg-2
                    checkj = checkj && !isassigned(poly, idx * babydeg + j)
                end
                checkj && continue

                # Compute the monomial if needed in the future.
                basis[(1<<i)+j] = mul(basis[j], basis[1<<i], rlk, oper)
                encode_to!(two, 2, oper, level=basis[(1<<i)+j].level[], scaling_factor=1.0)
                mul_to!(basis[(1<<i)+j], basis[(1<<i)+j], two, oper, islazy=true)

                # Set to the right precision.
                nowlvl, tarlvl = basis[(1<<i)-j].level[], basis[(1<<i)+j].level[]
                now_Qlen, now_auxQ = oper.Qatlevel[nowlvl+1]
                tar_Qlen, tar_auxQ = oper.Qatlevel[tarlvl+1]

                tmp = oper.ct_buff[3][1:now_Qlen]
                copy!(tmp, basis[(1<<i)-j])
                now_evalQ = geteval_at(now_Qlen, operQ, auxQ=now_auxQ)
                nowQ = prod(now_evalQ.moduli)

                @views tarQ = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
                Δ_real = basis[(1<<i)+j].scale[] * (nowQ // tarQ) / basis[(1<<i)-j].scale[]
                Δ_int = round(BigInt, Δ_real)

                for i = 1:now_Qlen
                    Δ_Qi = (Δ_int % now_evalQ[i].Q.Q) % UInt64
                    mul_to!(tmp.val.b.coeffs[i], Δ_Qi, tmp.val.b.coeffs[i], now_evalQ[i])
                    mul_to!(tmp.val.a.coeffs[i], Δ_Qi, tmp.val.a.coeffs[i], now_evalQ[i])
                end

                scale_to!(tmp.val, tmp.val, tar_Qlen, operQ, auxQ=tar_auxQ)
                tmp.level[] = tarlvl
                tmp.scale[] *= Δ_int / (nowQ // tarQ)

                # Subtract two ciphertexts.
                sub_to!(basis[(1<<i)+j], basis[(1<<i)+j], tmp, oper)
            end
        end
    end

    return basis
end

function get_scale_babystep(targetlvl::Integer, Δ_target::BigFloat, Ti::CKKS, oper::CKKSOperator)::BigFloat
    # Get the current level.
    nowlvl, Δ_Ti = Ti.level[], Ti.scale[]
    now_Qlen, now_auxQ = oper.Qatlevel[nowlvl+1]
    nowm1_Qlen, nowm1_auxQ = oper.Qatlevel[nowlvl]
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Compute Qdiff.
    evalQ = oper.operQ.evalQ
    @views nowQ = now_auxQ == 0 ? prod(evalQ.moduli[1:now_Qlen]) : prod(evalQ.moduli[1:now_Qlen-1]) * now_auxQ
    @views nowm1Q = nowm1_auxQ == 0 ? prod(evalQ.moduli[1:nowm1_Qlen]) : prod(evalQ.moduli[1:nowm1_Qlen-1]) * nowm1_auxQ
    @views tarQ = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ

    # Compute and output the scaling factor.
    Δ_target * nowQ / (Δ_Ti * tarQ * round(nowm1Q // tarQ))
end

function get_scale_giantstep(Δ_target::BigFloat, hilvl::Integer, lolvl::Integer, Ti::CKKS, oper::CKKSOperator)::Tuple{BigFloat,BigFloat}
    # Get the current level.
    Tilvl, Δ_Ti = Ti.level[], Ti.scale[]
    mullvl = min(Tilvl, hilvl)
    targetlvl = min(mullvl - 1, lolvl)

    Ti_Qlen, Ti_auxQ = oper.Qatlevel[Tilvl+1]
    hi_Qlen, hi_auxQ = oper.Qatlevel[hilvl+1]
    mul_Qlen, mul_auxQ = oper.Qatlevel[mullvl+1]

    lo_Qlen, lo_auxQ = oper.Qatlevel[lolvl+1]
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Compute Qdiff.
    evalQ = oper.operQ.evalQ

    @views QTi = Ti_auxQ == 0 ? prod(evalQ.moduli[1:Ti_Qlen]) : prod(evalQ.moduli[1:Ti_Qlen-1]) * Ti_auxQ
    @views Qhi = hi_auxQ == 0 ? prod(evalQ.moduli[1:hi_Qlen]) : prod(evalQ.moduli[1:hi_Qlen-1]) * hi_auxQ
    @views Qmul = mul_auxQ == 0 ? prod(evalQ.moduli[1:mul_Qlen]) : prod(evalQ.moduli[1:mul_Qlen-1]) * mul_auxQ

    Qdiff_real_Ti = QTi // Qmul
    Qdiff_int_Ti = round(Qdiff_real_Ti)

    Qdiff_real_hi = Qhi // Qmul
    Qdiff_int_hi = round(Qdiff_real_hi)

    @views Qlo = lo_auxQ == 0 ? prod(evalQ.moduli[1:lo_Qlen]) : prod(evalQ.moduli[1:lo_Qlen-1]) * lo_auxQ
    @views Qtar = tar_auxQ == 0 ? prod(evalQ.moduli[1:tar_Qlen]) : prod(evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ

    Qdiff_real_lo = Qlo // Qtar
    Qdiff_int_lo = round(Qdiff_real_lo)

    # Get the approximate scaling factor.
    mulm1_Qlen, mulm1_auxQ = oper.Qatlevel[mullvl]
    @views Qmulm1 = mulm1_auxQ == 0 ? prod(evalQ.moduli[1:mulm1_Qlen]) : prod(evalQ.moduli[1:mulm1_Qlen-1]) * mulm1_auxQ

    # Compute and output the scaling factor.
    Δ_lo = Δ_target * Qdiff_real_lo / Qdiff_int_lo
    Δ_hi = Δ_target * (Qmul // Qmulm1) / Δ_Ti * (Qdiff_real_hi / Qdiff_int_hi) * (Qdiff_real_Ti / Qdiff_int_Ti)

    Δ_lo, Δ_hi
end

function evalrecurse(lo::Int64, hi::Int64, poly::PolyCoeffs, basis::Dict{Int64,CKKS}, Δ_target::BigFloat, rlk::RLEV, oper::CKKSOperator)::Union{CKKS,Nothing}
    # Hyperparameters.
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    maxdeg = (1 << maxlevel) - 1
    babydeg = 1 << ceil(Int64, maxlevel / 2)

    # Current degree.
    d = hi - lo + 1
    if !ispow2(d)
        throw(DomainError("The degree should be a power of 2."))
    end

    # Corner cases.
    if lo > deg
        return nothing
    elseif d == 1
        res = similar(basis[1])
        if isassigned(poly, lo)
            coeffi = encode(poly[lo], oper, level=basis[1].level[] + 1, scaling_factor=Δ_target)
            add_to!(res, res, coeffi, oper)
            return res
        else
            return nothing
        end
    end

    # Temporary variables.
    ctlen = length(basis[1])
    tmp = oper.ct_buff[3][1:ctlen]

    if (d ≤ babydeg && hi < maxdeg) || (d == 2 && hi == maxdeg)
        coeffi = PlainConst(ctlen)

        # Compute the target level.
        targetlvl = basis[1].level[] - trailing_zeros(d) - 1
        hi == maxdeg && (targetlvl += 1)

        # Define result.
        N = oper.operQ.param.N
        tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]
        res = CKKS(N, tar_Qlen, targetlvl, auxQ=tar_auxQ, scale=Δ_target)

        # Babystep computation.
        if isassigned(poly, lo)
            encode_to!(coeffi, poly[lo], oper, level=res.level[], scaling_factor=Δ_target)
            add_to!(res, res, coeffi, oper)
        end

        for i = 1:min(d - 1, deg - lo)
            if isassigned(poly, lo + i)
                scale = get_scale_babystep(targetlvl, Δ_target, basis[i], oper)
                encode_to!(coeffi, poly[lo+i], oper, level=basis[i].level[], scaling_factor=scale)
                mul_to!(tmp, coeffi, basis[i], oper)
                add_to!(res, res, tmp, oper)
            end
        end

        return res
    else
        # Giantstep computation.
        halfd = d >> 1
        loghalfd = trailing_zeros(halfd)

        hilvl = basis[1].level[] - (hi == maxdeg ? loghalfd : loghalfd + 1)
        lolvl = basis[1].level[] - loghalfd - 1
        Δ_lo, Δ_hi = get_scale_giantstep(Δ_target, hilvl, lolvl, basis[halfd], oper)

        lo_res = evalrecurse(lo, lo + halfd - 1, poly, basis, Δ_lo, rlk, oper)
        hi_res = evalrecurse(lo + halfd, hi, poly, basis, Δ_hi, rlk, oper)
        isnothing(lo_res) && return nothing
        drop_level_to!(lo_res, lo_res, min(hilvl - 1, basis[halfd].level[] - 1, lolvl), oper)

        if !isnothing(hi_res)
            mul_to!(tmp, hi_res, basis[halfd], rlk, oper)
            add_to!(lo_res, lo_res, tmp, oper)
        end

        return lo_res
    end
end