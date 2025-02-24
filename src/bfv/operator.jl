"""
BFVOperator is a struct for arithmetic operations over BFV ciphertexts.
"""
struct BFVOperator
    ptxt_modulus::Modulus
    operQ::Operator
    evalR::PolyEvaluator
    packer::Union{IntPacker,Missing}
    ct_buff::Vector{BFV}
    tensor_buff::Tensor
    Qatlevel::Vector{Int64}
    Ratlevel::Vector{Int64}

    function BFVOperator(param::BFVParameters)
        ring_param, P, Q, dlen, t, ispacking, islevelled = param.ring_param, param.P, param.Q, param.dlen, param.ptxt_modulus, param.ispacking, param.islevelled

        # define operator.
        rlwe_param = RLWEParameters(ring_param, P, Q, dlen)
        operQ = Operator(rlwe_param)

        # define packer
        packer = ispacking ? IntPacker(t, ring_param) : missing

        # Setting R.
        Rprimes = collect(find_prime(ring_param, 62, 2length(Q) + 1))
        filter!(x -> x ∉ Q, Rprimes)
        Ridx, logR, logQ = 1, 0.0, log2(operQ.evalQ.moduli)
        while true
            logR += log2(Rprimes[Ridx])
            if logR > logQ + log2(ring_param.m)
                break
            end
            Ridx += 1
        end
        R = Modulus.(Rprimes[1:Ridx])
        evalR = PolyEvaluator(ring_param, R)

        # define buff.
        Qlen, Rlen = length(Q), length(R)
        ct_buff = BFV[BFV(ring_param.N, Qlen + Rlen, 0) for _ = 1:4]
        tensor_buff = Tensor(ring_param.N, Qlen + Rlen)

        # Compute the length of modulus chain.
        sfac = ring_param.m * big(√12) * t
        minQ = Float64(log2(2 * t * 6sfac))

        if islevelled
            Qatlevel = Vector{Int64}(undef, 0)
            Qlen = length(Q)
            auxQ = UInt64(0)
            while true
                pushfirst!(Qatlevel, Qlen)

                nowQ = 0.0
                @inbounds for i = 1:Qlen-1
                    nowQ += log2(Q[i])
                end
                nowQ += auxQ == 0 ? log2(Q[Qlen]) : log2(auxQ)
                nowQ < minQ && break

                auxQ == 0 && (auxQ = Q[Qlen])
                if sfac > auxQ
                    Qlen -= 1
                    auxQ = round(UInt64, 1 / sfac * Q[Qlen] * auxQ)
                    auxQ == 1 && (auxQ = UInt64(0))
                else
                    auxQ = round(UInt64, 1 / sfac * auxQ)
                    if auxQ == 1
                        auxQ = UInt64(0)
                        Qlen -= 1
                    end
                end
            end
        else
            maxlevel = floor(Int64, log2(operQ.evalQ.moduli) / log2(sfac))
            Qatlevel = fill(length(Q), maxlevel)
        end

        Ratlevel = similar(Qatlevel)
        @views for j = 1:length(Qatlevel)
            if islevelled
                Rlen, logR, logQ = 1, 0.0, log2(operQ.evalQ.moduli[1:Qatlevel[j]])
                while true
                    logR += log2(Rprimes[Rlen])
                    if logR > logQ + log2(ring_param.m)
                        break
                    end
                    Rlen += 1
                end
            else
                Rlen = length(R)
            end

            Ratlevel[j] = Rlen
        end

        new(Modulus(t), operQ, evalR, packer, ct_buff, tensor_buff, Qatlevel, Ratlevel)
    end

    function BFVOperator(oper::BFVOperator, ptxt_modulus::UInt64, top_Qlen::Int64; islevelled::Bool=true, ispacking::Bool=false)
        # Define others.
        operQ, evalR, ct_buff, tensor_buff, ring_param = oper.operQ, oper.evalR, oper.ct_buff, oper.tensor_buff, oper.operQ.param
        packer = ispacking ? IntPacker(ptxt_modulus, ring_param) : missing

        # Compute the length of modulus chain.

        sfac = ring_param.m * big(√12) * ptxt_modulus
        minQ = Float64(log2(2 * ptxt_modulus * 6sfac))

        if islevelled
            Qatlevel = Vector{Int64}(undef, 0)
            Q = operQ.evalQ.moduli
            Qlen, auxQ = top_Qlen, UInt64(0) #TODO better noise estimation?
            while true
                pushfirst!(Qatlevel, Qlen)

                nowQ = 0.0
                @inbounds for i = 1:Qlen-1
                    nowQ += log2(Q[i].Q)
                end
                nowQ += auxQ == 0 ? log2(Q[Qlen].Q) : log2(auxQ)
                nowQ < minQ && break

                auxQ == 0 && (auxQ = Q[Qlen].Q)
                if sfac > auxQ
                    Qlen -= 1
                    auxQ = round(UInt64, 1 / sfac * Q[Qlen].Q * auxQ)
                    auxQ == 1 && (auxQ = UInt64(0))
                else
                    auxQ = round(UInt64, 1 / sfac * auxQ)
                    if auxQ == 1
                        auxQ = UInt64(0)
                        Qlen -= 1
                    end
                end
            end
        else
            maxlevel = floor(Int64, log2(operQ.evalQ.moduli) / log2(sfac))
            Qatlevel = fill(length(operQ.evalQ.moduli), maxlevel)
        end

        Ratlevel = similar(Qatlevel)
        @views for j = 1:length(Qatlevel)
            if islevelled
                Rlen, logR, logQ = 1, 0.0, log2(operQ.evalQ.moduli[1:Qatlevel[j]])
                while true
                    logR += log2(Rprimes[Rlen])
                    if logR > logQ + log2(ring_param.m)
                        break
                    end
                    Rlen += 1
                end
            else
                Rlen = length(R)
            end

            Ratlevel[j] = Rlen
        end

        new(Modulus(ptxt_modulus), operQ, evalR, packer, ct_buff, tensor_buff, Qatlevel, Ratlevel)
    end
end

#======================================================================================#
################################## ENCODE AND DECODE ###################################
#======================================================================================#

function encode(msg::UInt64, oper::BFVOperator; level::Integer=typemax(Int64), isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen = Qatlevel[level+1]

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

function encode_to!(res::PlainConst, msg::UInt64, oper::BFVOperator; level::Integer=typemax(Int64), isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    if isPQ
        @assert !ismissing(operQ.evalP) "The parameter does not support PQ."

        Qlen = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ, isPQ=true)
        Plen = length(operQ.evalP)

        resize!(res, Plen + Qlen)
        for i = 1:Plen+Qlen
            res.val.vals[i] = _Bred(msg, evalQ[i])
        end
        res.auxQ[] = 0
        res.isPQ[] = true
    else
        Qlen = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ)

        resize!(res, Qlen)
        for i = 1:Qlen
            res.val.vals[i] = _Bred(msg, evalQ[i])
        end
        res.auxQ[] = 0
        res.isPQ[] = false
    end
end

function encode(msg::Vector{UInt64}, oper::BFVOperator; level::Integer=typemax(Int64), ispacking::Bool=true, isPQ::Bool=false)
    @assert level ≥ 0 "The level should be greater than 0."

    Qatlevel = oper.Qatlevel
    if level == typemax(Int64)
        level = length(Qatlevel) - 1
    end

    Qlen = Qatlevel[level+1]

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

function encode_to!(res::PlainPoly, msg::Vector{UInt64}, oper::BFVOperator; level::Integer=typemax(Int64), ispacking::Bool=true, isPQ::Bool=false)
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

        Qlen = Qatlevel[level+1]
        Plen = length(operQ.evalP)
        evalPQ = geteval_at(Plen + Qlen, operQ, isPQ=true)

        resize!(res, Plen + Qlen)
        be = BasisExtender([t], evalPQ.moduli)
        basis_extend!(res.val.coeffs, buff, be)
        res.val.isntt[] = false
        ntt!(res.val, evalPQ)

        res.auxQ[] = 0
        res.isPQ[] = true
    else
        Qlen = Qatlevel[level+1]
        evalQ = geteval_at(Qlen, operQ)

        resize!(res, Qlen)
        be = BasisExtender([t], evalQ.moduli)
        basis_extend!(res.val.coeffs, buff, be)
        res.val.isntt[] = false
        ntt!(res.val, evalQ)

        res.auxQ[] = 0
        res.isPQ[] = false
    end
end

function decode(x::PlainConst, oper::BFVOperator)
    isPQ = x.isPQ[]
    evalQ = geteval_at(length(x), oper.operQ, isPQ=isPQ)

    be = BasisExtender(evalQ.moduli, [oper.ptxt_modulus])
    @views res = oper.tensor_buff[1].coeffs[1][1:1]
    basis_extend!(res, x.val.vals, be)

    res[1]
end

function decode(x::PlainPoly, oper::BFVOperator; ispacking::Bool=true)
    if ispacking && !ismissing(oper.packer)
        res = Vector{UInt64}(undef, oper.packer.k)
    else
        res = Vector{UInt64}(undef, oper.operQ.param.N)
    end

    decode_to!(res, x, oper, ispacking=ispacking)

    res
end

function decode_to!(res::Vector{UInt64}, x::PlainPoly, oper::BFVOperator; ispacking::Bool=true)
    isPQ = x.isPQ[]
    evalQ = geteval_at(length(x), oper.operQ, isPQ=isPQ)

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

#======================================================================================#
################################## PLAINTEXT OPERATIONS ################################
#======================================================================================#

#TODO plaintext addition/subtraction/multiplication

change_level(x::PlainText, targetlvl::Integer, oper::BFVOperator) = begin
    res = similar(x)
    change_level_to!(res, x, targetlvl, oper)
    res
end

function change_level_to!(res::PlainPoly, x::PlainPoly, targetlvl::Integer, oper::BFVOperator)
    @assert targetlvl ≥ 0 "The target level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen = Qatlevel[targetlvl+1]

    change_modulus_to!(res, x, tar_Qlen, operQ)
end

function change_level_to!(res::PlainConst, x::PlainConst, targetlvl::Integer, oper::BFVOperator)
    @assert targetlvl ≥ 0 "The target level should be greater than 0."

    operQ, Qatlevel = oper.operQ, oper.Qatlevel
    tar_Qlen = Qatlevel[targetlvl+1]

    change_modulus_to!(res, x, tar_Qlen, operQ)
end

#================================================================================================#
################################## CIPHERTEXT OPERATIONS ########################################
#================================================================================================#

function drop_level_to!(res::BFV, x::BFV, targetlvl::Integer, oper::BFVOperator)
    currentlvl = x.level[]

    @assert targetlvl ≥ 0 "The target level should be greater than 0."
    @assert targetlvl ≤ currentlvl "The target level should be less than or equal to the current level."

    # define operator.
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    # Get the length of the modulus chain at the current level.
    now_Qlen, tar_Qlen = Qatlevel[currentlvl+1], Qatlevel[targetlvl+1]

    if now_Qlen == tar_Qlen
        resize!(res.val, length(x.val))
        copy!(res, x)
        return
    end

    # Sanity check
    @assert now_Qlen == length(x.val) "Something is wrong with the ciphertext."

    # Copy the ciphertext to the buffer, while dropping the unnecessary moduli for faster arithmetic.
    buff = oper.operQ.ct_buff[1][1:now_Qlen]
    copy!(buff, x.val)

    # Rational rescale to the new modulus.
    resize!(res.val, tar_Qlen)
    scale_to!(res.val, buff, tar_Qlen, operQ)

    # Update the level of the ciphertext.
    res.level[] = targetlvl
end

rescale(x::BFV, oper::BFVOperator) = begin
    res = similar(x)
    rescale_to!(res, x, oper)
    res
end

rescale_to!(res::BFV, x::BFV, oper::BFVOperator) = drop_level_to!(res, x, x.level[] - 1, oper)

function neg(x::BFV, oper::BFVOperator)
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::BFV, x::BFV, oper::BFVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    tar_Qlen = Qatlevel[level+1]
    evalQ = geteval_at(tar_Qlen, operQ)

    # Compute res = -x.
    resize!(res.val, tar_Qlen)
    for i = 1:tar_Qlen
        _neg_to!(res.val.b[i], x.val.b[i], evalQ[i])
        _neg_to!(res.val.a[i], x.val.a[i], evalQ[i])
    end
    res.level[] = level
end

add(x::BFV, y::PlainPoly, oper::BFVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainPoly, y::BFV, oper::BFVOperator) = add(y, x, oper)

function add_to!(res::BFV, x::BFV, y::PlainPoly, oper::BFVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == 0 "Something is wrong with the ciphertext."

    # Match the level.
    buff = PlainPoly(oper.tensor_buff[1][1:Qlen])
    change_level_to!(buff, y, level, oper)

    # Compute buff *= Δ
    evalQ = geteval_at(Qlen, operQ)
    Δ = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q)
    for i = eachindex(evalQ)
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        _mul_to!(buff.val[i], Δi, buff.val[i], evalQ[i])
    end

    # Compute x + buff.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
end

add_to!(res::BFV, x::PlainPoly, y::BFV, oper::BFVOperator) = add_to!(res, y, x, oper)

add(x::BFV, y::PlainConst, oper::BFVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add(x::PlainConst, y::BFV, oper::BFVOperator) = add(y, x, oper)

function add_to!(res::BFV, x::BFV, y::PlainConst, oper::BFVOperator)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == 0 "Something is wrong with the ciphertext."

    # Match the level.
    buff = change_level(y, level, oper)

    # Compute buff *= Δ
    evalQ = geteval_at(Qlen, operQ)
    Δ = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q)
    for i = eachindex(evalQ)
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        buff.val[i] = _mul(Δi, buff.val[i], evalQ[i])
    end

    # Compute x + m.
    resize!(res.val, Qlen)
    add_to!(res.val, x.val, buff, operQ)
end

add_to!(res::BFV, x::UInt64, y::BFV, oper::BFVOperator) = add_to!(res, y, x, oper)

add(x::BFV, y::BFV, oper::BFVOperator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

"""
Add two BFV ciphertexts.
"""
function add_to!(res::BFV, x::BFV, y::BFV, oper::BFVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)

    tar_Qlen = oper.Qatlevel[targetlvl+1]
    tmpx, tmpy = oper.ct_buff[1][1:tar_Qlen], oper.ct_buff[2][1:tar_Qlen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    resize!(res.val, tar_Qlen)
    add_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

sub(x::BFV, y::PlainPoly, oper::BFVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

@views function sub_to!(res::BFV, x::BFV, y::PlainPoly, oper::BFVOperator)
    buff = PlainPoly(oper.tensor_buff[2][1:length(x)])
    neg_to!(buff, y, oper.operQ)
    add_to!(res, x, buff, oper)
end

sub(x::PlainPoly, y::BFV, oper::BFVOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BFV, x::PlainPoly, y::BFV, oper::BFVOperator) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BFV, y::PlainConst, oper::BFVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::BFV, x::BFV, y::PlainConst, oper::BFVOperator)
    tmp = neg(y, oper.operQ)
    add_to!(res, x, tmp, oper)
end

sub(x::PlainConst, y::BFV, oper::BFVOperator) = begin
    res = similar(y)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::BFV, x::PlainConst, y::BFV, oper::BFVOperator) = begin
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

sub(x::BFV, y::BFV, oper::BFVOperator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

"""
Add sub BFV ciphertexts.
"""
function sub_to!(res::BFV, x::BFV, y::BFV, oper::BFVOperator)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    tar_Qlen = oper.Qatlevel[targetlvl+1]

    # Sanity check
    @assert length(res.val) == tar_Qlen "The length of the ciphertext should match the length of the modulus chain."

    # Match the levels.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.ct_buff[1][1:xlen], oper.ct_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    sub_to!(res.val, tmpx.val, tmpy.val, oper.operQ)
    res.level[] = targetlvl
end

mul(x::BFV, y::PlainPoly, oper::BFVOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainPoly, y::BFV, oper::BFVOperator; islazy::Bool=false) = mul(y, x, oper, islazy=islazy)

function mul_to!(res::BFV, x::BFV, y::PlainPoly, oper::BFVOperator; islazy::Bool=false)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == 0 "Something is wrong with the ciphertext."

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

mul_to!(res::BFV, x::PlainPoly, y::BFV, oper::BFVOperator; islazy::Bool=false) = mul_to!(res, y, x, oper, islazy=islazy)

mul(x::BFV, y::PlainConst, oper::BFVOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, oper, islazy=islazy)
    res
end

mul(x::PlainConst, y::BFV, oper::BFVOperator; islazy::Bool=false) = mul(y, x, oper, islazy=islazy)

function mul_to!(res::BFV, x::BFV, y::PlainConst, oper::BFVOperator; islazy::Bool=false)
    operQ, Qatlevel = oper.operQ, oper.Qatlevel

    level = x.level[]
    Qlen = Qatlevel[level+1]

    @assert length(x) == Qlen && x.val.auxQ[] == 0 "Something is wrong with the ciphertext."

    # Match the level.
    buff = change_level(y, level, oper)

    # Compute x * y.
    resize!(res.val, Qlen)
    mul_to!(res.val, x.val, buff, operQ)
    res.level[] = level

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)
end

mul_to!(res::BFV, x::PlainConst, y::BFV, oper::BFVOperator; islazy::Bool=false) = mul_to!(res, y, x, oper, islazy=islazy)

mul(x::BFV, y::BFV, rlk::RLEV, oper::BFVOperator; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, rlk, oper, islazy=islazy)
    res
end

function mul_to!(res::BFV, x::BFV, y::BFV, rlk::RLEV, oper::BFVOperator; islazy::Bool=false)
    xlevel, ylevel = x.level[], y.level[]
    targetlvl = min(xlevel, ylevel)
    @assert targetlvl > 0 "The input ciphertexts should be at least at level 1."

    # Drop the unnecessary levels of input ciphertexts.
    xlen, ylen = length(x.val), length(y.val)
    tmpx, tmpy = oper.ct_buff[1][1:xlen], oper.ct_buff[2][1:ylen]

    drop_level_to!(tmpx, x, targetlvl, oper)
    drop_level_to!(tmpy, y, targetlvl, oper)

    # Lift.
    Qlen, Rlen = oper.Qatlevel[targetlvl+1], oper.Ratlevel[targetlvl+1]
    liftlen = Qlen + Rlen
    liftx, lifty = oper.ct_buff[1][1:liftlen], oper.ct_buff[2][1:liftlen]

    _lift_to!(liftx, tmpx, oper)
    _lift_to!(lifty, tmpy, oper)

    # Tensoring.
    buffQR = oper.tensor_buff[1:3, 1:liftlen]
    _tensor_to!(buffQR, liftx, lifty, oper)

    # Division by Δ.
    evalQ, evalR = oper.operQ.evalQ[1:Qlen], oper.evalR[1:Rlen]
    evalQR = vcat(evalQ, evalR)
    intt_to!(buffQR.vals[1], buffQR.vals[1], evalQR)
    intt_to!(buffQR.vals[2], buffQR.vals[2], evalQR)
    intt_to!(buffQR.vals[3], buffQR.vals[3], evalQR)

    Q, R = evalQ.moduli, evalR.moduli
    cs = ComplexScaler(vcat(Q, R), Q, oper.ptxt_modulus.Q // prod(Q))

    buffQ = buffQR[1:3, 1:Qlen]
    complex_scale!(buffQ.vals[1].coeffs, buffQR.vals[1].coeffs, cs)
    complex_scale!(buffQ.vals[2].coeffs, buffQR.vals[2].coeffs, cs)
    complex_scale!(buffQ.vals[3].coeffs, buffQR.vals[3].coeffs, cs)

    # Relinearisation.
    resize!(res.val, Qlen)
    relinearise_to!(res.val, buffQ, rlk, oper.operQ)
    res.level[] = targetlvl

    # Rescale the ciphertext.
    !islazy && rescale_to!(res, res, oper)
end

function _lift_to!(res::BFV, x::BFV, oper::BFVOperator)
    xlevel = x.level[]

    # define operator.
    operQ = oper.operQ
    Qlen, Rlen = oper.Qatlevel[xlevel+1], oper.Ratlevel[xlevel+1]

    # Sanity check
    @assert length(x.val) == Qlen && length(res.val) == Qlen + Rlen "The length of the ciphertext should match the length of the modulus chain."

    # define evaluators.
    evalQ = operQ.evalQ[1:Qlen]
    evalR = oper.evalR[1:Rlen]

    # define BasisExtender.
    be = BasisExtender(evalQ.moduli, evalR.moduli)

    # Define buffer.
    buffQ = oper.ct_buff[3].val[1:Qlen]
    buffR = oper.ct_buff[3].val[Qlen+1:Qlen+Rlen]

    # Lifting.
    if x.val.b.isntt[]
        for i = 1:Qlen
            _intt_to!(buffQ.b.coeffs[i], x.val.b.coeffs[i], evalQ[i])
            _intt_to!(buffQ.a.coeffs[i], x.val.a.coeffs[i], evalQ[i])
            copy!(res.val.b.coeffs[i], x.val.b.coeffs[i])
            copy!(res.val.a.coeffs[i], x.val.a.coeffs[i])
        end
    else
        for i = 1:Qlen
            copy!(buffQ.b.coeffs[i], x.val.b.coeffs[i])
            copy!(buffQ.a.coeffs[i], x.val.a.coeffs[i])
            _ntt_to!(res.val.b.coeffs[i], x.val.b.coeffs[i], evalQ[i])
            _ntt_to!(res.val.a.coeffs[i], x.val.a.coeffs[i], evalQ[i])
        end
    end

    basis_extend!(buffR.b.coeffs, buffQ.b.coeffs, be)
    basis_extend!(buffR.a.coeffs, buffQ.a.coeffs, be)

    # Copy the result to the output in NTT form.
    for i = 1:Rlen
        _ntt_to!(res.val.b.coeffs[Qlen+i], buffR.b.coeffs[i], evalR[i])
        _ntt_to!(res.val.a.coeffs[Qlen+i], buffR.a.coeffs[i], evalR[i])
    end

    res.level[] = xlevel
    res.val.auxQ[] = 0
    res.val.isPQ[] = false
    res.val.b.isntt[] = true
    res.val.a.isntt[] = true
end

function _tensor_to!(res::Tensor, x::BFV, y::BFV, oper::BFVOperator)
    # Sanity check
    @assert x.level[] == y.level[] "The level of input ciphertexts should match."

    # define operator.
    lvl = x.level[]
    Qlen, Rlen = oper.Qatlevel[lvl+1], oper.Ratlevel[lvl+1]

    # Sanity check
    _, len = size(res)
    @assert len == length(x.val) == length(y.val) == Qlen + Rlen "The length of the ciphertext should match the parameters."

    # Tensor.
    evalQR = vcat(oper.operQ.evalQ[1:Qlen], oper.evalR[1:Rlen])

    mul_to!(res.vals[1], x.val.b, y.val.b, evalQR)
    mul_to!(res.vals[2], x.val.a, y.val.b, evalQR)
    muladd_to!(res.vals[2], x.val.b, y.val.a, evalQR)
    mul_to!(res.vals[3], x.val.a, y.val.a, evalQR)

    res.auxQ[] = 0
    res.isPQ[] = false
end

function decompose_a(x::BFV, oper::BFVOperator)
    len, decer = length(x.val), oper.operQ.decer

    if ismissing(oper.operQ.evalP)
        dlen = ceil(Int64, len / decer.dlen)
        deca = Tensor(oper.operQ.param.N, len, dlen)
    else
        dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.operQ.evalP)
        deca = Tensor(oper.operQ.param.N, len+Plen, dlen)
    end

    decompose_a_to!(deca, x, oper)

    deca
end

function decompose_a_to!(deca::Tensor, x::BFV, oper::BFVOperator)
    targetlvl = x.level[]

    # Sanity check
    Qatlevel = oper.Qatlevel
    tar_Qlen = Qatlevel[targetlvl+1]
    @assert length(x.val) == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Decompose.
    poly = PlainPoly(x.val.a)
    decompose_to!(deca, poly, oper.operQ)
end

keyswitch(x::BFV, ksk::RLEV, oper::BFVOperator) = begin
    res = similar(x)
    keyswitch_to!(res, x, ksk, oper)
    res
end

function keyswitch_to!(res::BFV, ct::BFV, ksk::RLEV, oper::BFVOperator)
    targetlvl = ct.level[]

    # Sanity check
    tar_Qlen = oper.Qatlevel[targetlvl+1]
    @assert length(ct.val) == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Copy the ciphertext to the buffer.
    buff = oper.ct_buff[1][1:tar_Qlen]

    # Key switching.
    resize!(res.val, tar_Qlen)
    keyswitch_to!(res.val, buff.val, ksk, oper.operQ)

    res.level[] = ct.level[]
end

rotate(x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BFVOperator) where {N} = begin
    res = similar(x)
    rotate_to!(res, x, idx, rtk, oper)
    res
end

function rotate_to!(res::BFV, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BFVOperator) where {N}
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, x, autidx, rtk, oper)
end

hoisted_rotate(adec::Tensor, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BFVOperator) where {N} = begin
    res = similar(x)
    hoisted_rotate_to!(res, adec, x, idx, rtk, oper)
    res
end

function hoisted_rotate_to!(res::BFV, adec::Tensor, x::BFV, idx::NTuple{N,Int64}, rtk::RLEV, oper::BFVOperator) where {N}
    packer = oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, x, autidx, rtk, oper)
end

automorphism(x::BFV, idx::Int64, rtk::RLEV, oper::BFVOperator) = begin
    res = similar(x)
    automorphism_to!(res, x, idx, rtk, oper)
    res
end

function automorphism_to!(res::BFV, x::BFV, idx::Int64, atk::RLEV, oper::BFVOperator)
    @assert length(res.val) == length(x.val) "The length of input and output ciphertexts should match."

    # Sanity check
    tar_Qlen = oper.Qatlevel[x.level[]+1]
    @assert length(x.val) == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Automorphism.
    resize!(res.val, tar_Qlen)
    automorphism_to!(res.val, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
end

hoisted_automorphism(adec::Tensor, x::BFV, idx::Int64, rtk::RLEV, oper::BFVOperator) = begin
    res = similar(x)
    hoisted_automorphism_to!(res, adec, x, idx, rtk, oper)
    res
end

function hoisted_automorphism_to!(res::BFV, adec::Tensor, x::BFV, idx::Int64, atk::RLEV, oper::BFVOperator)
    @assert length(res.val) == length(x.val) "The length of input and output ciphertexts should match."

    # Sanity check
    tar_Qlen = oper.Qatlevel[x.level[]+1]
    @assert length(x.val) == tar_Qlen "The length of the tensor should match the length of the ciphertexts."

    # Automorphism.
    resize!(res.val, tar_Qlen)
    hoisted_automorphism_to!(res.val, adec, x.val, idx, atk, oper.operQ)

    res.level[] = x.level[]
end

export BFVOperator