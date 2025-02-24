"""
CKKSOperator is a struct for arithmetic operations over CKKS ciphertexts.
"""
struct CKKSOperator
    scaling_factor::UInt128
    operQ::Operator
    packer::ComplexPacker
    pack_buff::Vector{Int128}
    ct_buff::Vector{CKKS}
    tensor_buff::Tensor
    Qatlevel::Vector{Tuple{Int64,UInt64}}

    function CKKSOperator(param::CKKSParameters)
        ring_param, P, Q, dlen, Δ = param.ring_param, param.P, param.Q, param.dlen, param.scaling_factor

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
        # TODO Do we have parameters where the gap is smaller than 3 bits?
        minQ = 3 + log2(Δ)
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
            if Δ > auxQ
                Qlen -= 1
                auxQ = round(UInt64, 1 / Δ * Q[Qlen] * auxQ)
                auxQ == 1 && (auxQ = UInt64(0))
            else
                auxQ = round(UInt64, 1 / Δ * auxQ)
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end

        new(Δ, operQ, packer, pack_buff, ct_buff, tensor_buff, Qatlevel)
    end

    # Used for plaintext switching.
    function CKKSOperator(oper::CKKSOperator, scaling_factor::Integer, top_Qlen::Int64, top_auxQ::PlainConst)
        # Define others.
        operQ, pack_buff, ct_buff, tensor_buff = oper.operQ, oper.pack_buff, oper.ct_buff, oper.tensor_buff
        packer = oper.packer

        # Compute the length of modulus chain, and auxQ at each level.
        minQ = Float64(Δ * 2^3)
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
                Qlen -= 1
                auxQ = round(UInt64, 1 / Δ * Q[Qlen] * auxQ)
                auxQ == 1 && (auxQ = UInt64(0))
            else
                auxQ = round(UInt64, 1 / Δ * auxQ)
                if auxQ == 1
                    auxQ = UInt64(0)
                    Qlen -= 1
                end
            end
        end

        new(scaling_factor, operQ, packer, pack_buff, ct_buff, tensor_buff, Qatlevel)
    end
end

#======================================================================================#
################################## ENCODE AND DECODE ###################################
#======================================================================================#

function encode(msg::AbstractFloat, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false)
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

function encode_to!(res::PlainConst, msg::AbstractFloat, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    operQ, Qatlevel, Δ = oper.operQ, oper.Qatlevel, oper.scaling_factor
    tmp = round(Int128, msg * Δ)

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
            res.val.vals[i] = _Bred(tmp, evalQ[i])
        end
        res.auxQ[] = auxQ
        res.isPQ[] = true
    else
        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)

        resize!(res, Qlen)
        for i = 1:Qlen
            res.val.vals[i] = _Bred(tmp, evalQ[i])
        end
        res.auxQ[] = auxQ
        res.isPQ[] = false
    end
end

function encode(msg::Vector{<:Union{AbstractFloat,Complex{<:AbstractFloat}}}, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true)
    @assert level ≥ 0 "The level should be greater than 0."

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen, _ = Qatlevel[level+1]

    if isPQ
        @assert !ismissing(oper.operQ.evalP) "The parameter does not support PQ."
        res = PlainPoly(oper.operQ.param.N, Qlen + length(oper.operQ.evalP))
        encode_to!(res, msg, oper, level=level, isPQ=true, ispacking=ispacking)
    else
        res = PlainPoly(oper.operQ.param.N, Qlen)
        encode_to!(res, msg, oper, level=level, isPQ=false, ispacking=ispacking)
    end

    res
end

# Coefficient Packing
function encode_to!(res::PlainPoly, msg::Vector{<:Union{AbstractFloat,Complex{<:AbstractFloat}}}, oper::CKKSOperator; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true)
    @assert level ≥ 0 "The level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    if ispacking
        packer, pack_buff, Δ = oper.packer, oper.pack_buff, oper.scaling_factor
        pack_to!(pack_buff, msg, Δ, packer)
    else
        @assert length(pack_buff) == length(msg) "The length of the message should match the ring size."
        @. pack_buff = round(Int128, msg * Δ)
    end

    if isPQ
        @assert !ismissing(operQ.evalP) "The parameter does not support PQ."

        Qlen, auxQ = Qatlevel[level+1]
        Plen = length(operQ.evalP)
        evalPQ = geteval_at(Plen + Qlen, operQ, auxQ=auxQ, isPQ=true)

        resize!(res, Plen + Qlen)
        for i = eachindex(evalPQ)
            _Bred_to!(res.val.coeffs[i], pack_buff, evalPQ[i].Q)
        end

        res.val.isntt[] = false
        ntt!(res.val, evalPQ)

        res.auxQ[] = auxQ
        res.isPQ[] = true
    else
        Qlen, auxQ = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, auxQ=auxQ)

        resize!(res, Qlen)
        for i = eachindex(evalQ)
            _Bred_to!(res.val.coeffs[i], pack_buff, evalQ[i].Q)
        end

        res.val.isntt[] = false
        ntt!(res.val, evalQ)

        res.auxQ[] = auxQ
        res.isPQ[] = false
    end
end

#TODO decode to Int128. Can it be done without BigInteger arithmetic?
function decode(x::PlainConst, oper::CKKSOperator)
    Double64(to_big(x, oper.operQ)) / oper.scaling_factor
end

function decode(x::PlainPoly, oper::CKKSOperator; ispacking::Bool=true)
    if ispacking
        res = similar(oper.packer.buff)
    else
        res = Vector{Double64}(undef, x.val.N)
    end

    decode_to!(res, x, oper, ispacking=ispacking)

    res
end

function decode_to!(res::Vector{<:Union{Double64,ComplexDF64}}, x::PlainPoly, oper::CKKSOperator; ispacking::Bool=true)
    isPQ, auxQ = x.isPQ[], x.auxQ[]
    evalQ = geteval_at(length(x), oper.operQ, auxQ=auxQ, isPQ=isPQ)

    buff = oper.tensor_buff.vals[1][1:length(x)]
    copy!(buff, x.val)
    buff.isntt[] && intt!(buff, evalQ)

    packer, pack_buff, Δ = oper.packer, oper.pack_buff, oper.scaling_factor
    pack_buff .= to_big(x, oper.operQ)

    if ispacking
        @assert length(res) == oper.packer.k "The length of the plaintext should match the ring size."
        unpack_to!(res, pack_buff, Δ, packer)
    else
        @assert length(res) == x.val.N "The length of the plaintext should match the ring size."
        @. res = pack_buff / Δ
    end
end

#=================================================================================================#
#################################### PLAINTEXT OPERATIONS #########################################
#=================================================================================================#

#TODO plaintext addtion/subtraction/multiplication

change_level(x::PlainText, targetlvl::Integer, oper::CKKSOperator) = begin
    res = similar(x)
    change_level_to!(res, x, targetlvl, oper)
    res
end

function change_level_to!(res::PlainPoly, x::PlainPoly, targetlvl::Integer, oper::CKKSOperator)
    @assert targetlvl ≥ 0 "The target level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    change_modulus_to!(res, x, tar_Qlen, operQ, auxQ=tar_auxQ)
end

function change_level_to!(res::PlainConst, x::PlainConst, targetlvl::Integer, oper::CKKSOperator)
    @assert targetlvl ≥ 0 "The target level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    change_modulus_to!(res, x, tar_Qlen, operQ, auxQ=tar_auxQ)
end

#================================================================================================#
#################################### CIPHERTEXT OPERATIONS ########################################
#================================================================================================#

drop_level(x::CKKS, targetlvl::Integer, oper::CKKSOperator) = begin
    res = similar(x)
    drop_level_to!(res, x, targetlvl, oper)
    res
end

#TODO
function drop_level_to!(res::CKKS, x::CKKS, targetlvl::Integer, oper::CKKSOperator)
    currentlvl = x.level[]

    @assert targetlvl ≥ 0 "The target level should be greater than 0."
    @assert targetlvl ≤ currentlvl "The target level should be less than or equal to the current level."

    if targetlvl == currentlvl
        resize!(res.val, length(x.val))
        copy!(res, x)
        return
    end

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert now_Qlen == length(x.val) && now_auxQ == x.val.auxQ[] "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer.
    buff = operQ.ct_buff[1][1:now_Qlen]
    copy!(buff, x.val[1:now_Qlen])

    # Multiply by the divided bits.
    evalQ = geteval_at(now_Qlen, operQ, auxQ=now_auxQ)
    @views tarQ = tar_auxQ == 0 ? prod(operQ.evalQ.moduli[1:tar_Qlen]) : prod(operQ.evalQ.moduli[1:tar_Qlen-1]) * tar_auxQ
    Qdiff = round(BigInt, prod(evalQ.moduli) // tarQ)
    for i = eachindex(evalQ)
        diffi = (Qdiff % evalQ[i].Q.Q) % UInt64
        _mul_to!(buff.b[i], diffi, buff.b[i], evalQ[i])
        _mul_to!(buff.a[i], diffi, buff.a[i], evalQ[i])
    end

    # Rational rescale to the new modulus.
    resize!(res.val, tar_Qlen)
    scale_to!(res.val, buff, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Update the level of the ciphertext.
    res.level[] = targetlvl
end

rescale(x::CKKS, oper::CKKSOperator) = begin
    res = similar(x)
    rescale_to!(res, x, oper)
    res
end

function rescale_to!(res::CKKS, x::CKKS, oper::CKKSOperator)
    currentlvl = x.level[]

    @assert currentlvl > 0 "The level of ciphertext should be greater than 0."

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus at the current level.
    now_Qlen, now_auxQ = Qatlevel[currentlvl+1]
    tar_Qlen, tar_auxQ = Qatlevel[currentlvl]

    # Sanity check
    @assert now_Qlen == length(x.val) && now_auxQ == x.val.auxQ[] "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer.
    buff = oper.ct_buff[1][1:now_Qlen]
    copy!(buff, x)

    # Rational rescale to the new modulus.
    scale_to!(res.val, buff.val, tar_Qlen, operQ, auxQ=tar_auxQ)

    # Update the level of the ciphertext.
    res.level[] = currentlvl - 1
end

function neg(x::CKKS, oper::CKKSOperator)
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::CKKS, x::CKKS, oper::CKKSOperator)
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

add(x::CKKS, y::PlainPoly, oper::CKKSOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainPoly, y::CKKS, oper::CKKSOperator) = add(y, x, oper)

function add_to!(res::CKKS, x::CKKS, y::PlainPoly, oper::CKKSOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == auxQ "Something is wrong with the ciphertext."

    # Match the level.
    buff = PlainPoly(oper.tensor_buff.vals[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute x + y.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
end

add_to!(res::CKKS, x::PlainPoly, y::CKKS, oper::CKKSOperator) = add_to!(res, y, x, oper)

add(x::CKKS, y::PlainConst, oper::CKKSOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainConst, y::CKKS, oper::CKKSOperator) = add(y, x, oper)

function add_to!(res::CKKS, x::CKKS, y::PlainConst, oper::CKKSOperator)
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

add_to!(res::CKKS, x::PlainConst, y::CKKS, oper::CKKSOperator) = add_to!(res, y, x, oper)

add(x::CKKS, y::CKKS, oper::CKKSOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

"""
Add two CKKS ciphertexts.
"""
function add_to!(res::CKKS, x::CKKS, y::CKKS, oper::CKKSOperator)
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

sub(x::CKKS, y::PlainPoly, oper::CKKSOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::CKKS, x::CKKS, y::PlainPoly, oper::CKKSOperator)
    buff = PlainPoly(oper.tensor_buff.vals[2][1:length(x)])
    neg_to!(buff, y, oper.operQ)
    add_to!(res, x, buff, oper)
end

sub(x::PlainPoly, y::CKKS, oper::CKKSOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::CKKS, x::PlainPoly, y::CKKS, oper::CKKSOperator)
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::CKKS, y::PlainConst, oper::CKKSOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::CKKS, x::CKKS, y::PlainConst, oper::CKKSOperator) = begin
    tmp = neg(y, oper.operQ)
    add_to!(res, x, tmp, oper)
end

sub(x::PlainConst, y::CKKS, oper::CKKSOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::CKKS, x::PlainConst, y::CKKS, oper::CKKSOperator) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::CKKS, y::CKKS, oper::CKKSOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

"""
Subtract two CKKS ciphertexts.
"""
function sub_to!(res::CKKS, x::CKKS, y::CKKS, oper::CKKSOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.ct_buff[1][1:xlen], oper.ct_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    sub_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

mul(x::CKKS, y::PlainPoly, oper::CKKSOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainPoly, y::CKKS, oper::CKKSOperator; islazy::Bool=false) = mul(y, x, oper, islazy=islazy)

function mul_to!(res::CKKS, x::CKKS, y::PlainPoly, oper::CKKSOperator; islazy::Bool=false)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen, auxQ = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == auxQ "Something is wrong with the ciphertext."

    # Match the level.
    buff = PlainPoly(oper.tensor_buff.vals[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute x * y.
    resize!(res.val, Qlen)
    mul_to!(res.val, x.val, buff, operQ)
    res.level[] = level

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)
end

mul_to!(res::CKKS, x::PlainPoly, y::CKKS, oper::CKKSOperator; islazy::Bool=false) = mul_to!(res, y, x, oper, islazy=islazy)

mul(x::CKKS, y::PlainConst, oper::CKKSOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainConst, y::CKKS, oper::CKKSOperator; islazy::Bool=false) = mul(y, x, oper, islazy=islazy)

function mul_to!(res::CKKS, x::CKKS, y::PlainConst, oper::CKKSOperator; islazy::Bool=false)
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

mul_to!(res::CKKS, x::PlainConst, y::CKKS, oper::CKKSOperator; islazy::Bool=false) = mul_to!(res, y, x, oper, islazy=islazy)

function mul(x::CKKS, y::CKKS, rlk::RLEV, oper::CKKSOperator; islazy::Bool=false)
    res = similar(x)
    mul_to!(res, x, y, rlk, oper, islazy=islazy)
    res
end

function mul_to!(res::CKKS, x::CKKS, y::CKKS, rlk::RLEV, oper::CKKSOperator; islazy::Bool=false)
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

function _tensor_to!(res::Tensor, x::CKKS, y::CKKS, oper::CKKSOperator)
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

function _relinearise_to!(res::RLWE, x::Tensor, rlk::RLEV, oper::CKKSOperator)
    # define operator.
    operQ = oper.operQ

    # Get the length of the modulus chain and the auxiliary modulus.
    _, Qlen = size(x)

    # Relinearise.
    resize!(res, Qlen)
    relinearise_to!(res, x, rlk, operQ)
end

function decompose_a(x::CKKS, oper::CKKSOperator)
    len, decer = length(x.val), oper.operQ.decer

    dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.operQ.evalP)
    deca = Tensor(oper.operQ.param.N, len + Plen, dlen)

    decompose_a_to!(deca, x, oper)

    deca
end

function decompose_a_to!(deca::Tensor, x::CKKS, oper::CKKSOperator)
    targetlvl = x.level[]

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."
    
    # Decompose.
    poly = PlainPoly(x.val.a, auxQ=tar_auxQ)

    decompose_to!(deca, poly, oper.operQ)
end

keyswitch(x::CKKS, ksk::RLEV, oper::CKKSOperator) = begin
    res = similar(x)
    keyswitch_to!(res, x, ksk, oper)
    res
end

function keyswitch_to!(res::CKKS, ct::CKKS, ksk::RLEV, oper::CKKSOperator)
    targetlvl = ct.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(ct.val) == tar_Qlen && ct.val.auxQ[] == tar_auxQ "Something is wrong with the ciphertext."

    # Key switching.
    resize!(res.val, tar_Qlen)
    keyswitch_to!(res.val, ct.val, ksk, oper.operQ)

    res.level[] = ct.level[]
end

rotate(x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, oper)
    res
end

function rotate_to!(res::CKKS, x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N}
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, oper)
end

hoisted_rotate(adec::Vector{ModPoly}, x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N} = begin
    res = similar(x)
    hoisted_rotate_to!(res, adec, x, idx, rtk, oper)
    res
end

hoisted_rotate_to!(res::CKKS, adec::Vector{ModPoly}, x::CKKS, idx::NTuple{N,Int64}, rtk::RLEV, oper::CKKSOperator) where {N} = begin
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, x, autidx, rtk, oper)
end

#TODO conjugation.

automorphism(x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator) = begin
    res = similar(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

function automorphism_to!(res::CKKS, x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator)
    targetlvl = x.level[]

    # Get the length of the modulus chain.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."

    # Automorphism.
    resize!(res.val, tar_Qlen)
    automorphism_to!(res.val, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
end

hoisted_automorphism(adec::Tensor, x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator) = begin
    res = similar(x)
    hoisted_automorphism_to!(res, adec, x, idx, atk, oper)
    res
end

function hoisted_automorphism_to!(res::CKKS, adec::Tensor, x::CKKS, idx::Int64, atk::RLEV, oper::CKKSOperator)
    targetlvl = x.level[]

    # Get the length of the modulus chain and the auxiliary modulus.
    tar_Qlen, tar_auxQ = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(x.val) == tar_Qlen && x.val.auxQ[] == tar_auxQ "The length of the tensor should match the length of the ciphertexts."

    # Automorphism.
    resize!(res.val, tar_Qlen)
    hoisted_automorphism_to!(res.val, adec, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
end

export CKKSOperator