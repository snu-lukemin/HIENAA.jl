"""
    BGVOperator(param::BGVParameters)

BGVOperator is a struct for arithmetic operations over BGV ciphertexts.
"""
struct BGVOperator
    ptxt_modulus::Modulus
    operQ::Operator
    packer::Union{IntPacker,Missing}
    ct_buff::Vector{BGV}
    tensor_buff::Tensor
    Qatlevel::Vector{Tuple{Int64,UInt64}}

    function BGVOperator(param::BGVParameters)
        ring_param, P, Q, dlen, t, ispacking = param.ring_param, param.P, param.Q, param.dlen, param.ptxt_modulus, param.ispacking

        # define operator.
        rlwe_param = RLWEParameters(ring_param, P, Q, dlen)
        operQ = Operator(rlwe_param)

        # define packer
        packer = ispacking ? IntPacker(t, ring_param) : missing

        # define buff.
        ct_buff = BGV[BGV(ring_param.N, length(Q), 0) for _ = 1:4]
        tensor_buff = Tensor(ring_param.N, length(Q) + length(P))

        # Compute the length of modulus chain, and auxQ at each level.
        sfac = ring_param.m * big(√12) * t
        minQ = Float64(log2(2 * t * 6sfac))
        Qatlevel = Vector{Tuple{Int64,UInt64}}(undef, 0)

        Qlen = length(Q)
        auxQ = UInt64(0)
        while true
            pushfirst!(Qatlevel, (Qlen, auxQ))

            nowQ = 0.0
            @inbounds for i = 1:Qlen-1
                nowQ += log2(Q[i])
            end
            nowQ += auxQ == 0 ? log2(Q[Qlen]) : log2(auxQ)
            nowQ < minQ && break

            auxQ == 0 && (auxQ = Q[Qlen])
            if sfac > auxQ
                Qlen -= 1
                auxQ = round(UInt64, 1 / (t * sfac) * Q[Qlen] * auxQ) * t + 1
                auxQ == 1 && (auxQ = UInt64(0))
            else
                auxQ = round(UInt64, 1 / (t * sfac) * auxQ) * t + 1
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end

        new(Modulus(t), operQ, packer, ct_buff, tensor_buff, Qatlevel)
    end

    # Used for plaintext switching.
    function BGVOperator(oper::BGVOperator, ptxt_modulus::Integer, top_Qlen::Int64, top_auxQ::PlainConst; ispacking::Bool=false)
        # Define others.
        operQ, ct_buff, tensor_buff, ring_param = oper.operQ, oper.ct_buff, oper.tensor_buff, oper.operQ.param
        packer = ispacking ? IntPacker(ptxt_modulus, ring_param) : missing

        # Compute the length of modulus chain, and auxQ at each level.
        sfac = ring_param.m * big(√12) * ptxt_modulus
        minQ = Float64(log2(2 * ptxt_modulus * 6sfac))
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
            if sfac > auxQ
                Qlen -= 1
                auxQ = round(UInt64, 1 / (ptxt_modulus * sfac) * Q[Qlen].Q * auxQ) * ptxt_modulus + 1
                auxQ == 1 && (auxQ = UInt64(0))
            else
                auxQ = round(UInt64, 1 / (ptxt_modulus * sfac) * auxQ) * ptxt_modulus + 1
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end

        new(Modulus(ptxt_modulus), operQ, packer, ct_buff, tensor_buff, Qatlevel)
    end
end

#======================================================================================#
################################## ENCODE AND DECODE ###################################
#======================================================================================#

function encode(msg::UInt64, oper::BGVOperator; level::Integer=typemax(Int64), isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen, _ = Qatlevel[level+1]

    if isPQ
        @assert !ismissing(oper.operQ.evalP) "The parameter does not support PQ."
        res = PlainConst(Qlen + length(oper.operQ.evalP))
        encode_to!(res, msg, oper, level=level, isPQ=true)
    else
        res = PlainConst(Qlen)
        encode_to!(res, msg, oper, level=level, isPQ=false)
    end

    res
end

function encode_to!(res::PlainConst, msg::UInt64, oper::BGVOperator; level::Integer=typemax(Int64), isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    if isPQ
        @assert !ismissing(operQ.evalP) "The parameter does not support PQ."

        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ, isPQ=true)
        Plen = length(operQ.evalP)

        resize!(res, Plen + Qlen)
        for i = 1:Plen+Qlen
            res.val.vals[i] = _Bred(msg, evalQ[i])
        end
        res.auxQ[] = auxQ
        res.isPQ[] = true
    else
        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)

        resize!(res, Qlen)
        for i = 1:Qlen
            res.val.vals[i] = _Bred(msg, evalQ[i])
        end
        res.auxQ[] = auxQ
        res.isPQ[] = false
    end
end

function encode(msg::Vector{UInt64}, oper::BGVOperator; level::Integer=typemax(Int64), ispacking::Bool=true, isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen, _ = Qatlevel[level+1]

    if isPQ
        @assert !ismissing(oper.operQ.evalP) "The parameter does not support PQ."
        res = PlainPoly(oper.operQ.param.N, Qlen + length(oper.operQ.evalP))
        encode_to!(res, msg, oper, level=level, ispacking=ispacking, isPQ=true)
    else
        res = PlainPoly(oper.operQ.param.N, Qlen)
        encode_to!(res, msg, oper, level=level, ispacking=ispacking, isPQ=false)
    end

    res
end

function encode_to!(res::PlainPoly, msg::Vector{UInt64}, oper::BGVOperator; level::Integer=typemax(Int64), ispacking::Bool=true, isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    t, packer, buff = oper.ptxt_modulus, oper.packer, oper.tensor_buff[1].coeffs[1]
    if !ismissing(packer) && ispacking
        pack_to!(buff, msg, packer)
    else
        @assert length(msg) == length(buff) "The length of the plaintext should match the ring size."
        @. buff = msg
    end

    buff = oper.tensor_buff[1].coeffs[1:1]

    if isPQ
        @assert !ismissing(operQ.evalP) "The parameter does not support PQ."

        Qlen, auxQ = Qatlevel[level+1]
        Plen = length(operQ.evalP)
        evalPQ = geteval_at(Plen + Qlen, operQ, auxQ=auxQ, isPQ=true)

        resize!(res, Plen + Qlen)
        be = BasisExtender([t], evalPQ.moduli)
        basis_extend!(res.val.coeffs, buff, be)
        res.val.isntt[] = false
        ntt!(res.val, evalPQ)

        res.auxQ[] = auxQ
        res.isPQ[] = true
    else
        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)

        resize!(res, Qlen)
        be = BasisExtender([t], evalQ.moduli)
        basis_extend!(res.val.coeffs, buff, be)
        res.val.isntt[] = false
        ntt!(res.val, evalQ)

        res.auxQ[] = auxQ
        res.isPQ[] = false
    end
end

function decode(x::PlainConst, oper::BGVOperator)
    isPQ, auxQ = x.isPQ[], x.auxQ[]
    evalQ = geteval_at(length(x), oper.operQ, auxQ=auxQ, isPQ=isPQ)

    be = BasisExtender(evalQ.moduli, [oper.ptxt_modulus])
    @views res = oper.tensor_buff[1].coeffs[1][1:1]
    basis_extend!(res, x.val.vals, be)

    res[1]
end

function decode(x::PlainPoly, oper::BGVOperator; ispacking::Bool=true)
    if ispacking && !ismissing(oper.packer)
        res = Vector{UInt64}(undef, oper.packer.k)
    else
        res = Vector{UInt64}(undef, oper.operQ.param.N)
    end

    decode_to!(res, x, oper, ispacking=ispacking)

    res
end

function decode_to!(res::Vector{UInt64}, x::PlainPoly, oper::BGVOperator; ispacking::Bool=true)
    isPQ, auxQ = x.isPQ[], x.auxQ[]
    evalQ = geteval_at(length(x), oper.operQ, auxQ=auxQ, isPQ=isPQ)

    buff = oper.tensor_buff[1][1:length(x)]
    copy!(buff, x.val)
    buff.isntt[] && intt!(buff, evalQ)

    be = BasisExtender(evalQ.moduli, [oper.ptxt_modulus])
    @views basis_extend!(buff.coeffs[1:1], buff.coeffs, be)

    if ispacking && !ismissing(oper.packer)
        @assert length(res) == oper.packer.k "The length of the plaintext should match the ring size."
        unpack_to!(res, buff.coeffs[1], oper.packer)
    else
        @assert length(res) == oper.operQ.param.N "The length of the plaintext should match the ring size."
        @. res = buff.coeffs[1]
    end
end

#=================================================================================================#
#################################### PLAINTEXT OPERATIONS #########################################
#=================================================================================================#

#TODO plaintext addtion/subtraction/multiplication

change_level(x::PlainText, targetlvl::Integer, oper::BGVOperator) = begin
    res = similar(x)
    change_level_to!(res, x, targetlvl, oper)
    res
end

function change_level_to!(res::PlainPoly, x::PlainPoly, targetlvl::Integer, oper::BGVOperator)
    @assert targetlvl ≥ 0 "The target level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    change_modulus_to!(res, x, tar_Qlen, operQ, auxQ=tar_auxQ)
end

function change_level_to!(res::PlainConst, x::PlainConst, targetlvl::Integer, oper::BGVOperator)
    @assert targetlvl ≥ 0 "The target level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    change_modulus_to!(res, x, tar_Qlen, operQ, auxQ=tar_auxQ)
end

#================================================================================================#
#################################### CIPHERTEXT OPERATIONS ########################################
#================================================================================================#

drop_level(x::BGV, targetlvl::Integer, oper::BGVOperator) = begin
    res = similar(x)
    drop_level_to!(res, x, targetlvl, oper)
    res
end

function drop_level_to!(res::BGV, x::BGV, targetlvl::Integer, oper::BGVOperator)
    currentlvl = x.level[]

    @assert targetlvl ≥ 0 "The target level should be greater than 0."
    @assert targetlvl ≤ currentlvl "The target level should be less than or equal to the current level."

    if targetlvl == currentlvl
        resize!(res.val, length(x.val))
        copy!(res, x)
        return
    end

    # define operator.
    t, operQ, Qatlevel = oper.ptxt_modulus.Q, oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert now_Qlen == length(x.val) && now_auxQ == x.val.auxQ[] "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer, while dropping the unnecessary moduli for faster arithmetic.
    buff = oper.operQ.ct_buff[1][1:tar_Qlen]
    copy!(buff, x.val[1:tar_Qlen])

    # Convert the ciphertext into a BFV ciphertext.
    # More precisely, we compute buff = (1-Q)/t * buff = [t⁻¹]_Q * buff.
    evalQ = tar_Qlen == now_Qlen ? geteval_at(tar_Qlen, operQ, auxQ=now_auxQ) : geteval_at(tar_Qlen, operQ)
    for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(buff.b[i], tinvQi, buff.b[i], evalQ[i])
        _mul_to!(buff.a[i], tinvQi, buff.a[i], evalQ[i])
    end

    # Rational rescale to the new modulus.
    resize!(res.val, tar_Qlen)
    scale_to!(res.val, buff, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Convert the ciphertext into BGV again.
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    for i = eachindex(evalQ)
        _mul_to!(res.val.b[i], t, res.val.b[i], evalQ[i])
        _mul_to!(res.val.a[i], t, res.val.a[i], evalQ[i])
    end

    # Update the level of the ciphertext.
    res.level[] = targetlvl
end

rescale(x::BGV, oper::BGVOperator) = begin
    res = similar(x)
    rescale_to!(res, x, oper)
    res
end

function rescale_to!(res::BGV, x::BGV, oper::BGVOperator)
    currentlvl = x.level[]

    @assert currentlvl > 0 "The level of ciphertext should be greater than 0."

    # define operator.
    t, operQ, Qatlevel = oper.ptxt_modulus.Q, oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[currentlvl]

    # Sanity check
    @assert now_Qlen == length(x.val) && now_auxQ == x.val.auxQ[] "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer.
    buff = oper.ct_buff[1][1:now_Qlen]
    copy!(buff, x)

    # Convert the ciphertext into a BFV ciphertext.
    # More precisely, we compute buff = (1-Q)/t * buff = [t⁻¹]_Q * buff.
    now_evalQ = geteval_at(now_Qlen, operQ, auxQ=now_auxQ)
    tar_evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    for i = eachindex(now_evalQ)
        tinvQi = invmod(t, now_evalQ[i].Q.Q)
        _mul_to!(buff.val.b[i], tinvQi, buff.val.b[i], now_evalQ[i])
        _mul_to!(buff.val.a[i], tinvQi, buff.val.a[i], now_evalQ[i])
    end

    # Rational rescale to the new modulus.
    scale_to!(res.val, buff.val, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Convert the ciphertext into BGV again.
    for i = eachindex(tar_evalQ)
        _mul_to!(res.val.b[i], t, res.val.b[i], tar_evalQ[i])
        _mul_to!(res.val.a[i], t, res.val.a[i], tar_evalQ[i])
    end

    # Update the level of the ciphertext.
    res.level[] = currentlvl - 1
end

function neg(x::BGV, oper::BGVOperator)
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::BGV, x::BGV, oper::BGVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen, tar_auxQ = Qatlevel[level+1]
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)

    # Compute res = -x.
    resize!(res.val, tar_Qlen)
    for i = 1:tar_Qlen
        _neg_to!(res.val.b[i], x.val.b[i], evalQ[i])
        _neg_to!(res.val.a[i], x.val.a[i], evalQ[i])
    end
    res.level[] = level
end

add(x::BGV, y::PlainPoly, oper::BGVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainPoly, y::BGV, oper::BGVOperator) = add(y, x, oper)

function add_to!(res::BGV, x::BGV, y::PlainPoly, oper::BGVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == auxQ "Something is wrong with the ciphertext."

    # Match the level.
    buff = PlainPoly(oper.tensor_buff[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute x + y.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
end

add_to!(res::BGV, x::PlainPoly, y::BGV, oper::BGVOperator) = add_to!(res, y, x, oper)

add(x::BGV, y::PlainConst, oper::BGVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainConst, y::BGV, oper::BGVOperator) = add(y, x, oper)

function add_to!(res::BGV, x::BGV, y::PlainConst, oper::BGVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]
    
    @assert length(x) == Qlen && x.val.auxQ[] == auxQ "Something is wrong with the ciphertext."

    # Match the level.
    buff = change_level(y, level, oper)

    # Compute x + y.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
end

add_to!(res::BGV, x::PlainConst, y::BGV, oper::BGVOperator) = add_to!(res, y, x, oper)

add(x::BGV, y::BGV, oper::BGVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

"""
Add two BGV ciphertexts.
"""
function add_to!(res::BGV, x::BGV, y::BGV, oper::BGVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    tar_Qlen, _ = oper.Qatlevel[targetlvl+1]
    tmpx, tmpy = oper.ct_buff[1][1:tar_Qlen], oper.ct_buff[2][1:tar_Qlen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    resize!(res.val, tar_Qlen)
    add_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

sub(x::BGV, y::PlainPoly, oper::BGVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::BGV, x::BGV, y::PlainPoly, oper::BGVOperator)
    buff = PlainPoly(oper.tensor_buff[2][1:length(x)])
    neg_to!(buff, y, oper.operQ)
    add_to!(res, x, buff, oper)
end

sub(x::PlainPoly, y::BGV, oper::BGVOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::BGV, x::PlainPoly, y::BGV, oper::BGVOperator)
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BGV, y::PlainConst, oper::BGVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BGV, x::BGV, y::PlainConst, oper::BGVOperator) = begin
    tmp = neg(y, oper.operQ)
    add_to!(res, x, tmp, oper)
end

sub(x::PlainConst, y::BGV, oper::BGVOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BGV, x::PlainConst, y::BGV, oper::BGVOperator) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BGV, y::BGV, oper::BGVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

"""
Subtract two BGV ciphertexts.
"""
function sub_to!(res::BGV, x::BGV, y::BGV, oper::BGVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.ct_buff[1][1:xlen], oper.ct_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    sub_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

mul(x::BGV, y::PlainPoly, oper::BGVOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainPoly, y::BGV, oper::BGVOperator; islazy::Bool=false) = mul(y, x, oper, islazy=islazy)

function mul_to!(res::BGV, x::BGV, y::PlainPoly, oper::BGVOperator; islazy::Bool=false)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]
    
    @assert length(x) == Qlen && x.val.auxQ[] == auxQ "Something is wrong with the ciphertext."

    # Match the level.
    buff = PlainPoly(oper.tensor_buff[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute x * y.
    resize!(res.val, Qlen)
    mul_to!(res.val, x.val, buff, operQ)
    res.level[] = level

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)
end

mul_to!(res::BGV, x::PlainPoly, y::BGV, oper::BGVOperator; islazy::Bool=false) = mul_to!(res, y, x, oper, islazy=islazy)

mul(x::BGV, y::PlainConst, oper::BGVOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainConst, y::BGV, oper::BGVOperator; islazy::Bool=false) = mul(y, x, oper, islazy=islazy)

function mul_to!(res::BGV, x::BGV, y::PlainConst, oper::BGVOperator; islazy::Bool=false)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]
    
    @assert length(x) == Qlen && x.val.auxQ[] == auxQ "Something is wrong with the ciphertext."

    # Match the level.
    buff = change_level(y, level, oper)

    # Compute x * y.
    resize!(res.val, Qlen)
    mul_to!(res.val, x.val, buff, operQ)
    res.level[] = level

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)
end

mul_to!(res::BGV, x::PlainConst, y::BGV, oper::BGVOperator; islazy::Bool=false) = mul_to!(res, y, x, oper, islazy=islazy)

function mul(x::BGV, y::BGV, rlk::RLEV, oper::BGVOperator; islazy::Bool=false)
    res = similar(x)
    mul_to!(res, x, y, rlk, oper, islazy=islazy)
    res
end

function mul_to!(res::BGV, x::BGV, y::BGV, rlk::RLEV, oper::BGVOperator; islazy::Bool=false)
    tar_Qlen = min(length(x.val), length(y.val))

    # Tensor the input ciphertexts.
    buff = oper.tensor_buff[1:3, 1:tar_Qlen]
    _tensor_to!(buff, x, y, oper)

    # Relinearisation.
    resize!(res.val, tar_Qlen)
    _relinearise_to!(res.val, buff, rlk, oper)
    res.level[] = min(x.level[], y.level[])

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)
end

function _tensor_to!(res::Tensor, x::BGV, y::BGV, oper::BGVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    @assert targetlvl > 0 "The input ciphertexts should be at least at level 1."

    # define operator
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the target level.
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    _, len = size(res)
    @assert len == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Drop the unnecessary levels of input ciphertexts.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.ct_buff[1][1:xlen], oper.ct_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    # Tensor.
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    !tmpx.val.b.isntt[] && ntt!(tmpx.val.b, evalQ)
    !tmpx.val.a.isntt[] && ntt!(tmpx.val.a, evalQ)
    !tmpy.val.b.isntt[] && ntt!(tmpy.val.b, evalQ)
    !tmpy.val.a.isntt[] && ntt!(tmpy.val.a, evalQ)

    mul_to!(res.vals[1], tmpx.val.b, tmpy.val.b, evalQ)
    mul_to!(res.vals[2], tmpx.val.a, tmpy.val.b, evalQ)
    muladd_to!(res.vals[2], tmpx.val.b, tmpy.val.a, evalQ)
    mul_to!(res.vals[3], tmpx.val.a, tmpy.val.a, evalQ)

    res.auxQ[] = tar_auxQ
end

function _relinearise_to!(res::RLWE, x::Tensor, rlk::RLEV, oper::BGVOperator)
    # define operator.
    t, operQ = oper.ptxt_modulus.Q, oper.operQ

    # Get the length of the modulus chain and the auxiliary modulus.
    _, Qlen = size(x)
    auxQ = x.auxQ[]

    # Copy the ciphertext to the buffer.
    buff = oper.tensor_buff[1:3, 1:Qlen]
    copy!(buff, x)

    # Convert the ciphertext into a BFV ciphertext.
    evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)
    for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(buff.vals[1][i], tinvQi, buff.vals[1][i], evalQ[i])
        _mul_to!(buff.vals[2][i], tinvQi, buff.vals[2][i], evalQ[i])
        _mul_to!(buff.vals[3][i], tinvQi, buff.vals[3][i], evalQ[i])
    end

    # Relinearise.
    resize!(res, Qlen)
    relinearise_to!(res, buff, rlk, operQ)

    # Convert back to BGV format.
    for i = eachindex(evalQ)
        _mul_to!(res.b.coeffs[i], t, res.b.coeffs[i], evalQ[i])
        _mul_to!(res.a.coeffs[i], t, res.a.coeffs[i], evalQ[i])
    end
end

function decompose_a(x::BGV, oper::BGVOperator)
    len, decer = length(x.val), oper.operQ.decer

    dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.operQ.evalP)
    deca = Tensor(oper.operQ.param.N, len+Plen, dlen)

    decompose_a_to!(deca, x, oper)

    deca
end

function decompose_a_to!(deca::Tensor, x::BGV, oper::BGVOperator)
    targetlvl = x.level[]

    # define operator.
    t, operQ, Qatlevel = oper.ptxt_modulus.Q, oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."

    # Copy the ciphertext to the buffer.
    buff = PlainPoly(oper.tensor_buff[1][1:tar_Qlen], auxQ=tar_auxQ)
    copy!(buff.val, x.val.a)

    # Convert the ciphertext into a BFV ciphertext.
    evalQ = geteval_at(tar_Qlen, operQ, auxQ=tar_auxQ)
    for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(buff.val[i], tinvQi, buff.val[i], evalQ[i])
    end

    decompose_to!(deca, buff, oper.operQ)
end

function _MSB_packing_to!(res::RLWE, ct::RLWE, oper::BGVOperator)
    # define operator.
    t, operQ = oper.ptxt_modulus.Q, oper.operQ

    # Sanity check.
    Qlen, auxQ = length(ct), ct.auxQ[]
    @assert length(res) == Qlen "The length of the tensor should match the length of the ciphertexts."
    @assert !ct.isPQ[] "The ciphertext should not be in PQ."

    # Convert the ciphertext into a BFV ciphertext.
    evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)
    for i = eachindex(evalQ)
        tinvQi = invmod(t, evalQ[i].Q.Q)
        _mul_to!(res.b[i], tinvQi, ct.b[i], evalQ[i])
        _mul_to!(res.a[i], tinvQi, ct.a[i], evalQ[i])
    end

    res.auxQ[] = auxQ
    res.isPQ[] = false
end

function _LSB_packing_to!(res::RLWE, ct::RLWE, oper::BGVOperator)
    # define operator.
    t, operQ = oper.ptxt_modulus.Q, oper.operQ

    # Sanity check.
    Qlen, auxQ = length(ct), ct.auxQ[]
    @assert length(res) == Qlen "The length of the tensor should match the length of the ciphertexts."
    @assert !ct.isPQ[] "The ciphertext should not be in PQ."

    # Convert the ciphertext into a BFV ciphertext.
    evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)
    for i = eachindex(evalQ)
        _mul_to!(res.b[i], t, ct.b[i], evalQ[i])
        _mul_to!(res.a[i], t, ct.a[i], evalQ[i])
    end

    res.auxQ[] = auxQ
    res.isPQ[] = false
end

keyswitch(x::BGV, ksk::RLEV, oper::BGVOperator) = begin
    res = similar(x)
    keyswitch_to!(res, x, ksk, oper)
    res
end

function keyswitch_to!(res::BGV, ct::BGV, ksk::RLEV, oper::BGVOperator)
    targetlvl = ct.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(ct.val) == tar_Qlen && ct.val.auxQ[] == tar_auxQ "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer.
    buff = oper.ct_buff[1][1:tar_Qlen]

    # Convert the ciphertext into a BFV ciphertext.
    _MSB_packing_to!(buff.val, ct.val, oper)

    # Key switching.
    resize!(res.val, tar_Qlen)
    keyswitch_to!(res.val, buff.val, ksk, oper.operQ)

    # Convert back to BGV format.
    _LSB_packing_to!(res.val, res.val, oper)

    res.level[] = ct.level[]
end

rotate(x::BGV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BGVOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, oper)
    res
end

function rotate_to!(res::BGV, x::BGV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BGVOperator) where {N}
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, oper)
end

hoisted_rotate(adec::Vector{ModPoly}, x::BGV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BGVOperator) where {N} = begin
    res = similar(x)
    hoisted_rotate_to!(res, adec, x, idx, rtk, oper)
    res
end

hoisted_rotate_to!(res::BGV, adec::Vector{ModPoly}, x::BGV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BGVOperator) where {N} = begin
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, x, autidx, rtk, oper)
end

automorphism(x::BGV, idx::Int64, atk::RLEV, oper::BGVOperator) = begin
    res = similar(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

function automorphism_to!(res::BGV, x::BGV, idx::Int64, atk::RLEV, oper::BGVOperator)
    targetlvl = x.level[]

    # Get the length of the modulus chain.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."

    # Convert the ciphertext into a BFV ciphertext.
    buff = oper.ct_buff[1][1:tar_Qlen]
    _MSB_packing_to!(buff.val, x.val, oper)

    # Automorphism.
    resize!(res.val, tar_Qlen)
    automorphism_to!(res.val, buff.val, idx, atk, oper.operQ)

    # Convert back to BGV format.
    _LSB_packing_to!(res.val, res.val, oper)

    res.level[] = x.level[]
end

hoisted_automorphism(adec::Tensor, x::BGV, idx::Int64, atk::RLEV, oper::BGVOperator) = begin
    res = similar(x)
    hoisted_automorphism_to!(res, adec, x, idx, atk, oper)
    res
end

function hoisted_automorphism_to!(res::BGV, adec::Tensor, x::BGV, idx::Int64, atk::RLEV, oper::BGVOperator)
    targetlvl = x.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."

    # Convert the ciphertext into a BFV ciphertext.
    buff = oper.ct_buff[1][1:tar_Qlen]
    _MSB_packing_to!(buff.val, x.val, oper)

    # Automorphism.
    resize!(res.val, tar_Qlen)
    hoisted_automorphism_to!(res.val, adec, buff.val, idx, atk, oper.operQ)

    # Convert back to BGV format.
    _LSB_packing_to!(res.val, res.val, oper)

    res.level[] = x.level[]
end

export BGVOperator