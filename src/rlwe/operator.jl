struct Operator
    param::RingParam
    evalP::Union{Missing,PolyEvaluator}
    evalQ::PolyEvaluator
    auxeval::_PolyEvaluatorWord
    decer::Decomposer
    ct_buff::Vector{RLWE}
    tensor_buff::Tensor

    function Operator(param::RLWEParameters)
        ring_param, P, Q, dlen = param.ring_param, param.P, param.Q, param.dlen
        N = ring_param.N

        if ismissing(P)
            Qlen = length(Q)
            Qmoduli = Modulus.(Q)
            evalP, evalQ = missing, PolyEvaluator(ring_param, Qmoduli)
            auxeval = _PolyEvaluatorWord(ring_param, Modulus(1 << 62))
            decer = Decomposer(Qmoduli, dlen)
            ct_buff = [RLWE(N, Qlen, isPQ=false) for _ = 1:2]
            tensor_buff = Tensor(N, Qlen, decer.glen + 2, isPQ=true)
        else
            Plen, Qlen = length(P), length(Q)
            Pmoduli, Qmoduli = Modulus.(P), Modulus.(Q)
            evalP, evalQ = PolyEvaluator(ring_param, Pmoduli), PolyEvaluator(ring_param, Qmoduli)
            auxeval = _PolyEvaluatorWord(ring_param, Modulus(1 << 62))
            decer = Decomposer(Pmoduli, Qmoduli, dlen)
            ct_buff = [RLWE(N, Plen + Qlen, isPQ=true) for _ = 1:2]
            tensor_buff = Tensor(N, Plen + Qlen, decer.glen + 2, isPQ=true)
        end

        new(ring_param, evalP, evalQ, auxeval, decer, ct_buff, tensor_buff)
    end
end

function geteval_at(len::Int64, oper::Operator; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
    evalP, evalQ = oper.evalP, oper.evalQ
    if isPQ
        @assert !ismissing(evalP) "Special modulus is not defined for the operator."
        Plen, Qlen = length(evalP), length(evalQ)
        @assert Qlen + Plen ≥ len > Plen
        eval = auxQ == 0 ? vcat(evalP, evalQ[1:len-Plen]) : vcat(evalP, evalQ[1:len-Plen-1], _PolyEvaluatorWord(oper.auxeval, Modulus(auxQ)))
    else
        Qlen = length(evalQ)
        @assert len ≤ Qlen
        eval = auxQ == 0 ? evalQ[1:len] : vcat(evalQ[1:len-1], _PolyEvaluatorWord(oper.auxeval, Modulus(auxQ)))
    end

    eval
end

#============================================================================================#
############################## Plaintext Operations ##########################################
#============================================================================================#

ntt(x::PlainPoly, oper::Operator) = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::PlainPoly, x::PlainPoly, oper::Operator) = begin
    @assert length(x) == length(res) "The input and output plaintexts have different length."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper; isPQ=isPQ, auxQ=auxQ)
    ntt_to!(res.val, x.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

intt(x::PlainPoly, oper::Operator) = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::PlainPoly, x::PlainPoly, oper::Operator) = begin
    @assert length(x) == length(res) "The input and output plaintexts have different length."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper; isPQ=isPQ, auxQ=auxQ)
    intt_to!(res.val, x.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

neg(x::PlainText, oper::Operator) = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

neg_to!(res::T, x::T, oper::Operator) where {T<:PlainText} = begin
    @assert length(res) == length(x) "The plaintexts should have the same length."

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    neg_to!(res.val, x.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
end

add_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator) = begin
    @assert length(res) == length(x) == length(y) "The plaintexts should have the same length."
    @assert x.isPQ[] == y.isPQ[] "The plaintexts should be all in PQ or not in PQ."
    @assert x.auxQ[] == y.auxQ[] "The plaintexts should have the same auxiliary modulus."

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    add_to!(res.val, x.val, y.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
end

sub_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator) = begin
    @assert length(res) == length(x) == length(y) "The plaintexts should have the same length."
    @assert x.isPQ[] == y.isPQ[] "The plaintexts should be all in PQ or not in PQ."
    @assert x.auxQ[] == y.auxQ[] "The plaintexts should have the same auxiliary modulus."

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    sub_to!(res.val, x.val, y.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
end

mul_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator) = begin
    @assert length(res) == length(x) == length(y) "The plaintexts should have the same length."
    @assert x.isPQ[] == y.isPQ[] "The plaintexts should be all in PQ or not in PQ."
    @assert x.auxQ[] == y.auxQ[] "The plaintexts should have the same auxiliary modulus."

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    mul_to!(res.val, x.val, y.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
end

muladd_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator) = begin
    @assert length(res) == length(x) == length(y) "The plaintexts should have the same length."
    @assert res.isPQ[] == x.isPQ[] == y.isPQ[] "The plaintexts should be all in PQ or not in PQ."
    @assert res.auxQ[] == x.auxQ[] == y.auxQ[] "The plaintexts should have the same auxiliary modulus."

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    muladd_to!(res.val, x.val, y.val, eval)
end

mulsub_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator) = begin
    @assert length(res) == length(x) == length(y) "The plaintexts should have the same length."
    @assert res.isPQ[] == x.isPQ[] == y.isPQ[] "The plaintexts should be all in PQ or not in PQ."
    @assert res.auxQ[] == x.auxQ[] == y.auxQ[] "The plaintexts should have the same auxiliary modulus."

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    mulsub_to!(res.val, x.val, y.val, eval)
end

to_big(x::PlainText, oper::Operator) = begin
    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    to_big(x.val, eval)
end

change_modulus(x::PlainText, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0)) = begin
    res = similar(x)
    change_modulus_to!(res, x, len, oper, auxQ=auxQ)
    res
end

function change_modulus_to!(res::PlainConst, x::PlainConst, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))
    isPQ = x.isPQ[]
    now_Qlen, now_auxQ = length(x), x.auxQ[]
    tar_Qlen, tar_auxQ = len, auxQ

    if isPQ
        @assert !ismissing(oper.evalP) "The parameter does not support PQ."
        tar_Qlen += length(oper.evalP)
    end

    if now_Qlen == tar_Qlen && now_auxQ == tar_auxQ
        resize!(res, length(x))
        copy!(res, x)
        return
    end

    eval_nowQ = geteval_at(now_Qlen, oper, auxQ=now_auxQ, isPQ=isPQ)
    eval_tarQ = geteval_at(tar_Qlen, oper, auxQ=tar_auxQ, isPQ=isPQ)

    # BasisExtend.
    be = BasisExtender(eval_nowQ.moduli, eval_tarQ.moduli)
    if tar_Qlen ≥ now_Qlen
        resize!(res, tar_Qlen)
        @views basis_extend!(res.val.vals[1:tar_Qlen], x.val.vals[1:now_Qlen], be)
    else
        @views basis_extend!(res.val.vals[1:tar_Qlen], x.val.vals[1:now_Qlen], be)
        resize!(res, tar_Qlen)
    end

    # Set parameters.
    res.auxQ[] = tar_auxQ
    res.isPQ[] = x.isPQ[]
end

function change_modulus_to!(res::PlainPoly, x::PlainPoly, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))
    isPQ = x.isPQ[]
    now_Qlen, now_auxQ = length(x), x.auxQ[]
    tar_Qlen, tar_auxQ = len, auxQ

    if isPQ
        @assert !ismissing(oper.evalP) "The parameter does not support PQ."
        tar_Qlen += length(oper.evalP)
    end

    if now_Qlen ≥ tar_Qlen && now_auxQ == tar_auxQ
        if now_auxQ == 0
            for i = 1:tar_Qlen
                copy!(res.val.coeffs[i], x.val.coeffs[i])
            end
        else
            for i = 1:tar_Qlen-1
                copy!(res.val.coeffs[i], x.val.coeffs[i])
            end
            copy!(res.val.coeffs[tar_Qlen], x.val.coeffs[end])
        end

        resize!(res, tar_Qlen)
        res.auxQ[] = tar_auxQ
        res.isPQ[] = x.isPQ[]
        res.val.isntt[] = x.val.isntt[]

        return
    end

    eval_nowQ = geteval_at(now_Qlen, oper, auxQ=now_auxQ, isPQ=isPQ)
    eval_tarQ = geteval_at(tar_Qlen, oper, auxQ=tar_auxQ, isPQ=isPQ)

    #TODO optimisation.

    # Copy the ciphertext to the buffer.
    buff = oper.tensor_buff[1][1:now_Qlen]
    buff2 = oper.tensor_buff[2][1:now_Qlen]
    copy!(buff, x.val)
    copy!(buff2, x.val)
    resize!(res, tar_Qlen)

    # BasisExtend.
    buff.isntt[] && intt!(buff, eval_nowQ)
    be = BasisExtender(eval_nowQ.moduli, eval_tarQ.moduli)
    basis_extend!(res.val.coeffs, buff.coeffs, be)

    # NTT, with optimisation.
    if buff2.isntt[]
        reuselen = min(now_Qlen, tar_Qlen) - 1
        for i = 1:reuselen
            copy!(res.val.coeffs[i], buff2.coeffs[i])
        end
        for i = reuselen+1:tar_Qlen
            _ntt!(res.val.coeffs[i], eval_tarQ[i])
        end
    end
    res.val.isntt[] = buff2.isntt[]

    # Set parameters.
    res.auxQ[] = tar_auxQ
    res.isPQ[] = x.isPQ[]
end

scale(x::PlainPoly, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0)) = begin
    res = similar(x)
    scale_to!(res, x, len, oper, auxQ=auxQ)
    res
end

function scale_to!(res::PlainPoly, x::PlainPoly, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    newlen, newauxQ = len, auxQ
    oldlen, oldauxQ = length(x), x.auxQ[]
    @assert newlen ≤ len "The new length should be less than or equal to the original length."

    if oldlen == newlen && oldauxQ == newauxQ
        copy!(res, x)
        return
    end

    # Define Moduli.
    Rlen = min(oldlen, newlen) - 1
    if Rlen > 0
        # RT -> RS
        evalRT = geteval_at(oldlen, oper, auxQ=oldauxQ)
        evalRS = geteval_at(newlen, oper, auxQ=newauxQ)

        evalR, evalT = evalRT[1:Rlen], evalRT[Rlen+1:end]
        buffRT = oper.tensor_buff[1, 1:oldlen]
        buffRS = oper.tensor_buff[1, 1:newlen]
        buffR, buffT = buffRT[1:Rlen], buffRT[Rlen+1:end]

        R, T, RS = evalR.moduli, evalT.moduli, evalRS.moduli
        ss = SimpleScaler(T, RS)

        # Compute RinvT and TinvS.
        Rinv = ModScalar(1, evalT)
        TinvS = ModScalar(1, evalR)

        for i = eachindex(R), j = eachindex(T)
            Rinv.vals[j] = _Bmul(Rinv.vals[j], invmod(R[i].Q, T[j].Q), T[j])
            TinvS.vals[i] = _Bmul(TinvS.vals[i], invmod(T[j].Q, R[i].Q), R[i])
        end

        for i = eachindex(R), j = Rlen+1:newlen
            TinvS.vals[i] = _Bmul(TinvS.vals[i], RS[j].Q, R[i])
        end

        # scale b.
        copy!(buffRT, x.val)
        buffR.isntt[] = x.val.isntt[]
        buffT.isntt[] = x.val.isntt[]
        buffRS.isntt[] = x.val.isntt[]
        resize!(res, newlen)

        mul_to!(buffR, TinvS, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.val.coeffs, buffT.coeffs, ss)
        res.val.isntt[] = false

        buffRS.isntt[] && ntt!(res.val, evalRS)
        @. buffRS.coeffs[end] = 0
        add_to!(res.val, res.val, buffRS, evalRS)
    else
        # Define Moduli.
        # R -> T
        evalR = geteval_at(oldlen, oper, auxQ=oldauxQ)
        evalT = geteval_at(newlen, oper, auxQ=newauxQ)

        ss = SimpleScaler(evalR.moduli, evalT.moduli)
        buff = oper.tensor_buff[1, 1:oldlen]

        copy!(buff, x.val)
        buff.isntt[] && intt!(buff, evalR)
        simple_scale!((@view res.val.coeffs[1:newlen]), buff.coeffs, ss)
        resize!(res, newlen)
        res.val.isntt[] = false
        buff.isntt[] && ntt!(res.val, evalT)
    end

    res.isPQ[] = false
    res.auxQ[] = newauxQ
end

divide_and_round(x::PlainPoly, Δ::Real, oper::Operator) = begin
    res = similar(x)
    divide_and_round_to!(res, x, Δ, oper)
    res
end

"""
    Compute res = ⌊x (mod Q) / Δ⌉ (mod ⌊Q / Δ⌉).
"""
function divide_and_round_to!(res::PlainPoly, x::PlainPoly, Δ::Real, oper::Operator)
    @assert Δ ≥ 1 "The scaling factor should be greater than or equal to 1."
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(x), x.auxQ[]

    # Set Q.
    if auxQ == 0
        Q = oper.evalQ.moduli[1:len]
        auxQ = Q[end].Q
    else
        Q = vcat(oper.evalQ.moduli[1:len-1], Modulus(auxQ))
    end

    # Find the next Modulus.
    newlen = len
    setprecision(BigFloat, 64 * len)
    Δ = BigFloat(Δ)
    while Δ > auxQ
        newlen -= 1
        Δ /= Q[newlen].Q
    end

    newauxQ = round(UInt64, auxQ / Δ)

    if newauxQ == 1
        newlen -= 1
        newauxQ = UInt64(0)
    end

    scale_to!(res, x, newlen, oper, auxQ=newauxQ)
end

#============================================================================================#
############################## RLWE Operations ###############################################
#============================================================================================#

ntt(x::RLWE, oper::Operator) = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

function ntt_to!(res::RLWE, x::RLWE, oper::Operator)
    @assert length(x) == length(res) "The input and output ciphertexts have different length."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    ntt_to!(res.b, x.b, eval)
    ntt_to!(res.a, x.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

intt(x::RLWE, oper::Operator) = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

function intt_to!(res::RLWE, x::RLWE, oper::Operator)
    @assert length(x) == length(res) "The input and output ciphertexts have different length."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    intt_to!(res.b, x.b, eval)
    intt_to!(res.a, x.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

neg(x::RLWE, oper::Operator) = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::RLWE, x::RLWE, oper::Operator)
    @assert length(x) == length(res) "The input and output ciphertexts have different length."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    neg_to!(res.b, x.b, eval)
    neg_to!(res.a, x.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

add(x::T, y::S, oper::Operator) where {T,S<:Union{RLWE,PlainText}} = begin
    res = T == RLWE ? similar(x) : similar(y)
    add_to!(res, x, y, oper)
    res
end

function add_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    add_to!(res.b, x.val, y.b, eval)
    copy!(res.a, y.a)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

add_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator) = add_to!(res, y, x, oper)

function add_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    add_to!(res.b, x.b, y.b, eval)
    add_to!(res.a, x.a, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

sub(x::T, y::S, oper::Operator) where {T,S<:Union{RLWE,PlainText}} = begin
    res = T == RLWE ? similar(x) : similar(y)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator)
    @assert length(res) == length(x) == length(y) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    sub_to!(res.b, x.b, y.val, eval)
    copy!(res.a, x.a)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

function sub_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
end

function sub_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    sub_to!(res.b, x.b, y.b, eval)
    sub_to!(res.a, x.a, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

mul(x::T, y::S, oper::Operator) where {T,S<:Union{RLWE,PlainText}} = begin
    res = T == RLWE ? similar(x) : similar(y)
    mul_to!(res, x, y, oper)
    res
end

function mul_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)
    @assert length(x) == length(res) == length(y) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    mul_to!(res.b, x.val, y.b, eval)
    mul_to!(res.a, x.val, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

mul_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator) = mul_to!(res, y, x, oper)

function muladd_to!(res::RLWE, x::PlainText, y::PlainText, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    muladd_to!(res.b, x.val, y.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

function muladd_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    muladd_to!(res.b, x.val, y.b, eval)
    muladd_to!(res.a, x.val, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

muladd_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator) = muladd_to!(res, y, x, oper)

function mulsub_to!(res::RLWE, x::PlainText, y::PlainText, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] == res.isPQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    mulsub_to!(res.b, x.val, y.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

function mulsub_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different length."
    @assert x.auxQ[] == y.auxQ[] == res.isPQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    mulsub_to!(res.b, x.val, y.b, eval)
    mulsub_to!(res.a, x.val, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

mulsub_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator) = mulsub_to!(res, y, x, oper)

divide_by_P(x::RLWE, oper::Operator) = begin
    res = similar(x)
    divide_by_P_to!(res, x, oper)
    res
end

function divide_by_P_to!(res::RLWE, x::RLWE, oper::Operator)
    @assert x.isPQ[] "The input ciphertext should be in PQ."

    len, Plen, auxQ = length(x) - length(oper.evalP), length(oper.evalP), x.auxQ[]

    # Define Moduli.
    evalP, evalQ = oper.evalP, geteval_at(len, oper, auxQ=auxQ)
    buff = oper.tensor_buff[1, 1:Plen+len]
    buffP, buffQ = buff[1:Plen], buff[Plen+1:end]

    P, Q = evalP.moduli, evalQ.moduli
    ss = SimpleScaler(P, Q)

    # Compute Pinv and Qinv.
    Pinv = ModScalar(1, evalQ)
    Qinv = ModScalar(1, evalP)

    for i = eachindex(P), j = eachindex(Q)
        Pinv.vals[j] = _Bmul(Pinv.vals[j], invmod(P[i].Q, Q[j].Q), Q[j])
        Qinv.vals[i] = _Bmul(Qinv.vals[i], invmod(Q[j].Q, P[i].Q), P[i])
    end

    # scale b.
    copy!(buff, x.b)
    buffP.isntt[] = x.b.isntt[]
    buffQ.isntt[] = x.b.isntt[]

    mul_to!(buffP, Qinv, buffP, evalP)
    mul_to!(buffQ, Pinv, buffQ, evalQ)
    buffP.isntt[] && intt!(buffP, evalP)

    resize!(res.b, len)
    simple_scale!(res.b.coeffs, buffP.coeffs, ss)
    res.b.isntt[] = false

    buffQ.isntt[] && ntt!(res.b, evalQ)
    add_to!(res.b, res.b, buffQ, evalQ)

    # scale a.
    copy!(buff, x.a)
    buffP.isntt[] = x.a.isntt[]
    buffQ.isntt[] = x.a.isntt[]

    mul_to!(buffP, Qinv, buffP, evalP)
    mul_to!(buffQ, Pinv, buffQ, evalQ)
    buffP.isntt[] && intt!(buffP, evalP)

    resize!(res.a, len)
    simple_scale!(res.a.coeffs, buffP.coeffs, ss)
    res.a.isntt[] = false

    buffQ.isntt[] && ntt!(res.a, evalQ)
    add_to!(res.a, res.a, buffQ, evalQ)

    res.isPQ[] = false
    res.auxQ[] = auxQ
end

scale(x::RLWE, len::Int64; oper::Operator, auxQ::UInt64=UInt64(0)) = begin
    res = similar(x)
    scale_to!(res, x, len, oper, auxQ=auxQ)
    res
end

function scale_to!(res::RLWE, x::RLWE, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    newlen, newauxQ = len, auxQ
    oldlen, oldauxQ = length(x), x.auxQ[]
    @assert newlen ≤ len "The new length should be less than or equal to the original length."

    if oldlen == newlen && oldauxQ == newauxQ
        copy!(res, x)
        return
    end

    # Define Moduli.
    Rlen = min(oldlen, newlen) - 1
    if Rlen > 0
        # RT -> RS
        evalRT = geteval_at(oldlen, oper, auxQ=oldauxQ)
        evalRS = geteval_at(newlen, oper, auxQ=newauxQ)

        evalR, evalT = evalRT[1:Rlen], evalRT[Rlen+1:end]
        buffRT = oper.tensor_buff[1, 1:oldlen]
        buffRS = oper.tensor_buff[1, 1:newlen]
        buffR, buffT = buffRT[1:Rlen], buffRT[Rlen+1:end]

        R, T, RS = evalR.moduli, evalT.moduli, evalRS.moduli
        ss = SimpleScaler(T, RS)

        # Compute RinvT and TinvS.
        Rinv = ModScalar(1, evalT)
        TinvS = ModScalar(1, evalR)

        for i = eachindex(R), j = eachindex(T)
            Rinv.vals[j] = _Bmul(Rinv.vals[j], invmod(R[i].Q, T[j].Q), T[j])
            TinvS.vals[i] = _Bmul(TinvS.vals[i], invmod(T[j].Q, R[i].Q), R[i])
        end

        for i = eachindex(R), j = Rlen+1:newlen
            TinvS.vals[i] = _Bmul(TinvS.vals[i], RS[j].Q, R[i])
        end

        # scale b.
        copy!(buffRT, x.b)
        buffR.isntt[] = x.b.isntt[]
        buffT.isntt[] = x.b.isntt[]
        buffRS.isntt[] = x.b.isntt[]
        resize!(res.b, newlen)

        mul_to!(buffR, TinvS, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.b.coeffs, buffT.coeffs, ss)
        res.b.isntt[] = false

        buffRS.isntt[] && ntt!(res.b, evalRS)
        @. buffRS.coeffs[end] = 0
        add_to!(res.b, res.b, buffRS, evalRS)

        # scale a.
        copy!(buffRT, x.a)
        buffR.isntt[] = x.a.isntt[]
        buffT.isntt[] = x.a.isntt[]
        buffRS.isntt[] = x.a.isntt[]
        resize!(res.a, newlen)

        mul_to!(buffR, TinvS, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale!(res.a.coeffs, buffT.coeffs, ss)
        res.a.isntt[] = false

        buffRS.isntt[] && ntt!(res.a, evalRS)
        @. buffRS.coeffs[end] = 0
        add_to!(res.a, res.a, buffRS, evalRS)
    else
        # Define Moduli.
        # R -> T
        evalR = geteval_at(oldlen, oper, auxQ=oldauxQ)
        evalT = geteval_at(newlen, oper, auxQ=newauxQ)

        ss = SimpleScaler(evalR.moduli, evalT.moduli)
        buff = oper.tensor_buff[1, 1:oldlen]

        copy!(buff, x.b)
        buff.isntt[] && intt!(buff, evalR)
        simple_scale!((@view res.b.coeffs[1:newlen]), buff.coeffs, ss)
        resize!(res.b, newlen)
        res.b.isntt[] = false
        buff.isntt[] && ntt!(res.b, evalT)

        copy!(buff, x.a)
        buff.isntt[] && intt!(buff, evalR)
        simple_scale!((@view res.a.coeffs[1:newlen]), buff.coeffs, ss)
        resize!(res.a, newlen)
        res.a.isntt[] = false
        buff.isntt[] && ntt!(res.a, evalT)
    end

    res.isPQ[] = false
    res.auxQ[] = newauxQ
end

divide_and_round(x::RLWE, Δ::Real, oper::Operator) = begin
    res = similar(x)
    divide_and_round_to!(res, x, Δ, oper)
    res
end

"""
    Compute res = ⌊x (mod Q) / Δ⌉ (mod ⌊Q / Δ⌉).
"""
function divide_and_round_to!(res::RLWE, x::RLWE, Δ::Real, oper::Operator)
    @assert Δ ≥ 1 "The scaling factor should be greater than or equal to 1."
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(x), x.auxQ[]

    # Set Q.
    if auxQ == 0
        Q = oper.evalQ.moduli[1:len]
        auxQ = Q[end].Q
    else
        Q = vcat(oper.evalQ.moduli[1:len-1], Modulus(auxQ))
    end

    # Find the next Modulus.
    newlen = len
    setprecision(BigFloat, 64 * len)
    Δ = BigFloat(Δ)
    while Δ > auxQ
        newlen -= 1
        Δ /= Q[newlen].Q
    end

    newauxQ = round(UInt64, auxQ / Δ)

    if newauxQ == 1
        newlen -= 1
        newauxQ = UInt64(0)
    end

    scale_to!(res, x, newlen, oper, auxQ=newauxQ)
end

#============================================================================================#
############################## Tensor Operations #############################################
#============================================================================================#

ntt(x::Tensor, oper::Operator) = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::Tensor, x::Tensor, oper::Operator) = begin
    @assert length(x) == length(res) "The input and output ciphertexts have different sizes."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        ntt_to!(res.val[i], x.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

intt(x::Tensor, oper::Operator) = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::Tensor, x::Tensor, oper::Operator) = begin
    @assert length(x) == length(res) "The input and output ciphertexts have different sizes."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        intt_to!(res.val[i], x.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

neg(x::Tensor, oper::Operator) = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::Tensor, x::Tensor, oper::Operator)
    @assert length(x) == length(res) "The input and output ciphertexts have different sizes."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        neg_to!(res.val[i], x.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

add(x::Tensor, y::Tensor, oper::Operator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

function add_to!(res::Tensor, x::Tensor, y::Tensor, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        add_to!(res.val[i], x.val[i], y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

sub(x::Tensor, y::Tensor, oper::Operator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::Tensor, x::Tensor, y::Tensor, oper::Operator)
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        sub_to!(res.val[i], x.val[i], y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

mul(x::Tensor, y::PlainText, oper::Operator) = begin
    res = similar(x)
    mul_to!(res, x, y, oper)
    res
end

mul(x::PlainText, y::Tensor, oper::Operator) = mul(y, x, oper)

function mul_to!(res::Tensor, x::PlainText, y::Tensor, oper::Operator)
    @assert length(x) == length(res) == length(y) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] && x.isPQ[] == y.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        mul_to!(res.val[i], x.val, y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

mul_to!(res::Tensor, x::Tensor, y::PlainText, oper::Operator) = mul_to!(res, y, x, oper)

muladd_to!(res::Tensor, x::PlainText, y::Tensor, oper::Operator) = begin
    @assert length(x) == length(res) == length(y) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        muladd_to!(res.val[i], x.val, y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

muladd_to!(res::Tensor, x::Tensor, y::PlainText, oper::Operator) = muladd_to!(res, y, x, oper)

mulsub_to!(res::Tensor, x::PlainText, y::Tensor, oper::Operator) = begin
    @assert length(x) == length(y) == length(res) "The input and output ciphertexts have different sizes."
    @assert x.auxQ[] == y.auxQ[] == res.isPQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[] "The input and output ciphertexts have different moduli."

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        mulsub_to!(res.val[i], x.val, y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
end

mulsub_to!(res::Tensor, x::Tensor, y::PlainText, oper::Operator) = mulsub_to!(res, y, x, oper)

#================================================================================================#
##################################### RLEV Operations ############################################
#================================================================================================#

ntt(x::RLEV, oper::Operator) = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::RLEV, x::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen "The length of input and output ciphertext should match."

    for i = 1:res.len
        ntt_to!(res.stack[i], x.stack[i], oper)
    end
end

intt(x::RLEV, oper::Operator) = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::RLEV, x::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen "The length of input and output ciphertext should match."

    for i = 1:res.len
        intt_to!(res.stack[i], x.stack[i], oper)
    end
end

neg(x::RLEV, oper::Operator) = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

neg_to!(res::RLEV, x::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen "The length of input and output ciphertext should match."

    for i = 1:res.len
        neg_to!(res.stack[i], x.stack[i], oper)
    end
end

add(x::RLEV, y::RLEV, oper::Operator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add_to!(res::RLEV, x::RLEV, y::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen == y.glen "The length of input and output ciphertext should match."

    for i = 1:res.len
        add_to!(res.stack[i], x.stack[i], y.stack[i], oper.eval)
    end
end

sub(x::RLEV, y::RLEV, oper::Operator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::RLEV, x::RLEV, y::RLEV, oper::Operator) = begin
    @assert res.glen == x.glen == y.glen "The length of input and output ciphertext should match."

    for i = 1:res.len
        sub_to!(res.stack[i], x.stack[i], y.stack[i], oper.eval)
    end
end

#================================================================================================#
##################################### RGSW Operations ############################################
#================================================================================================#

ntt(x::RGSW, oper::Operator) = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::RGSW, x::RGSW, oper::Operator) = begin
    ntt_to!(res.basketb, x.basketb, oper)
    ntt_to!(res.basketa, x.basketa, oper)
end

intt(x::RGSW, oper::Operator) = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::RGSW, x::RGSW, oper::Operator) = begin
    intt_to!(res.basketb, x.basketb, oper)
    intt_to!(res.basketa, x.basketa, oper)
end

add(x::RGSW, y::RGSW, oper::Operator) = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add_to!(res::RGSW, x::RGSW, y::RGSW, oper::Operator) = begin
    add_to!(res.basketb, x.basketb, y.basketb, oper.eval)
    add_to!(res.basketa, x.basketa, y.basketa, oper.eval)
end

sub(x::RGSW, y::RGSW, oper::Operator) = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::RGSW, x::RGSW, y::RGSW, oper::Operator) = begin
    sub_to!(res.basketb, x.basketb, y.basketb, oper.eval)
    sub_to!(res.basketa, x.basketa, y.basketa, oper.eval)
end

#================================================================================================#
##################################### Decompositions #############################################
#================================================================================================#

decompose(x::PlainPoly, oper::Operator) = begin
    len, decer = length(x), oper.decer

    if ismissing(oper.evalP)
        dlen = ceil(Int64, len / decer.dlen)
        decx = Tensor(oper.param.N, len, dlen)
    else
        dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.evalP)
        decx = Tensor(oper.param.N, len + Plen, dlen)
    end

    decompose_to!(decx, x, oper)
    decx
end

function decompose_to!(decx::Tensor, x::PlainPoly, oper::Operator)
    @assert !x.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ, decer = length(x), x.auxQ[], oper.decer
    buff = PlainPoly(oper.tensor_buff[end, 1:len])

    copy!(buff, x)
    buff.val.isntt[] && intt_to!(buff, buff, oper)
    auxQ ≠ 0 && scale_to!(buff, buff, len, oper)

    buff.isPQ[] = false
    buff.auxQ[] = auxQ

    _decompose_to!(decx, buff, decer)
end

function hoisted_gadgetprod_to!(res::RLWE, decx::Tensor, ct::RLEV, oper::Operator; islazy::Bool=false)
    N, len = size(decx)
    !ismissing(oper.evalP) && (len -= length(oper.evalP))
    evalQ = geteval_at(len, oper)

    isPQ, auxQ = decx.isPQ[], decx.auxQ[]
    if !isPQ
        ct_buff = oper.ct_buff[1][1:len]

        initialise!(ct_buff, isntt=true, isPQ=false)

        for i = 1:N
            !decx[i].isntt[] && ntt_to!(decx[i], decx[i], evalQ)
            muladd_to!(ct_buff.b, decx[i], ct.stack[i].b[1:len], evalQ)
            muladd_to!(ct_buff.a, decx[i], ct.stack[i].a[1:len], evalQ)
        end

        resize!(res, len)
        copy!(res, ct_buff)
    else
        @assert !ismissing(oper.evalP) "The ciphertext is not defined over PQ."

        Plen = length(oper.evalP)
        evalPQ = geteval_at(len + Plen, oper, isPQ=true)
        ct_buff = oper.ct_buff[1][1:len+Plen]

        initialise!(ct_buff, isntt=true, isPQ=true)

        for i = 1:N
            !decx[i].isntt[] && ntt_to!(decx[i], decx[i], evalPQ)
            muladd_to!(ct_buff.b, decx[i], ct.stack[i].b[1:Plen+len], evalPQ)
            muladd_to!(ct_buff.a, decx[i], ct.stack[i].a[1:Plen+len], evalPQ)
        end

        if auxQ ≠ 0
            intt!(ct_buff.b, evalPQ)
            intt!(ct_buff.a, evalPQ)
        end

        if islazy
            resize!(res, length(ct_buff))
            copy!(res, ct_buff)

            return
        end

        divide_by_P_to!(res, ct_buff, oper)
    end

    if auxQ ≠ 0
        divide_and_round_to!(res, res, evalQ.moduli[end].Q // auxQ, oper)
        evalQ = geteval_at(len, oper, auxQ=auxQ)
    end

    !res.b.isntt[] && ntt!(res.b, evalQ)
    !res.a.isntt[] && ntt!(res.a, evalQ)
end

gadgetprod(x::PlainPoly, ct::RLEV, oper::Operator) = begin
    res = similar(ct.stack[1])
    gadgetprod_to!(res, x, ct, oper)
    res
end

function gadgetprod_to!(res::RLWE, x::PlainPoly, ct::RLEV, oper::Operator)
    len, decer = length(x), oper.decer

    if ismissing(oper.evalP)
        dlen = ceil(Int64, len / decer.dlen)
        decbuff = oper.tensor_buff[1:dlen, 1:len]
    else
        dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.evalP)
        decbuff = oper.tensor_buff[1:dlen, 1:len+Plen]
    end

    decompose_to!(decbuff, x, oper)
    hoisted_gadgetprod_to!(res, decbuff, ct, oper)
end

relinearise(ct::Tensor, rlk::RLEV, oper::Operator) = begin
    res = RLWE(similar(ct.vals[1]), similar(ct.vals[1]))
    relinearise_to!(res, ct, rlk, oper)
    res
end

function relinearise_to!(res::RLWE, ct::Tensor, rlk::RLEV, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."
    @assert length(ct.vals) == 3

    _, len = size(ct)
    auxQ = ct.auxQ[]
    evalQ = geteval_at(len, oper, auxQ=auxQ)
    buff = oper.tensor_buff[end, 1:len]

    gadgetprod_to!(res, PlainPoly(ct.vals[3], auxQ=auxQ), rlk, oper)

    copy!(buff, ct.vals[1])
    !buff.isntt[] && ntt!(buff, evalQ)
    add_to!(res.b, res.b, buff, evalQ)

    copy!(buff, ct.vals[2])
    !buff.isntt[] && ntt!(buff, evalQ)
    add_to!(res.a, res.a, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]
end

keyswitch(ct::RLWE, ksk::RLEV, oper::Operator) = begin
    res = similar(ct)
    keyswitch_to!(res, ct, ksk, oper)
    res
end

function keyswitch_to!(res::RLWE, ct::RLWE, ksk::RLEV, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = geteval_at(len, oper, auxQ=auxQ)
    buff = oper.tensor_buff[end-1, 1:len]

    copy!(buff, ct.b)
    !buff.isntt[] && ntt!(buff, evalQ)

    gadgetprod_to!(res, PlainPoly(ct.a, auxQ=auxQ), ksk, oper)

    add_to!(res.b, res.b, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]
end

hoisted_keyswitch(adec::Tensor, ct::RLWE, ksk::RLEV, oper::Operator) = begin
    res = similar(ct)
    hoisted_keyswitch_to!(res, adec, ct, ksk, oper)
    res
end

function hoisted_keyswitch_to!(res::RLWE, adec::Tensor, ct::RLWE, ksk::RLEV, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = geteval_at(len, oper, auxQ=auxQ)
    buff = oper.tensor_buff[end-1, 1:len]

    copy!(buff, ct.b)
    !buff.isntt[] && ntt!(buff, evalQ)

    hoisted_gadgetprod_to!(res, adec, ksk, oper)

    add_to!(res.b, res.b, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]
end

automorphism(x::RLWE, idx::Integer, atk::RLEV, oper::Operator) = begin
    res = deepcopy(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

automorphism!(x::RLWE, idx::Integer, atk::RLEV, oper::Operator) = automorphism_to!(x, x, idx, atk, oper)

function automorphism_to!(res::RLWE, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)
    eval = geteval_at(length(ct), oper, isPQ=ct.isPQ[], auxQ=ct.auxQ[])

    keyswitch_to!(res, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)
end

function hoisted_automorphism_to!(res::RLWE, adec::Tensor, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)
    eval = geteval_at(length(ct), oper, isPQ=ct.isPQ[], auxQ=ct.auxQ[])

    hoisted_keyswitch_to!(res, adec, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)
end

#======================================================================================#

function extprod(ct::RLWE, rgsw::RGSW, oper::Operator)
    res = similar(ct)
    extprod_to!(res, ct, rgsw, oper)
    res
end

function extprod_to!(res::RLWE, ct::RLWE, rgsw::RGSW, oper::Operator)
    @assert !ct.isPQ[] "The input ciphertext should not be in PQ."

    len, auxQ = length(ct), ct.auxQ[]
    buff = oper.tensor_buff[end-1, 1:len]
    ct_buff = oper.ct_buff[1][1:len]

    copy!(buff, ct.a)
    gadgetprod_to!(res, PlainPoly(ct.b, auxQ=auxQ), rgsw.basketb, oper)
    gadgetprod_to!(ct_buff, PlainPoly(buff, auxQ=auxQ), rgsw.basketa, oper)
    add_to!(res, res, ct_buff, oper)
end

@views function hoisted_extprod_to!(res::RLWE, ctdec::Tensor, rgsw::RGSW, oper::Operator)
    len = length(ctdec)
    !ismissing(oper.evalP) && (len -= length(oper.evalP))
    ct_buff, glen = oper.ct_buff[1][1:len], length(ctdec) >> 1

    hoisted_gadgetprod_to!(res, ctdec[1:glen], rgsw.basketb, oper)
    hoisted_gadgetprod_to!(ct_buff, ctdec[glen+1:end], rgsw.basketa, oper)
    add_to!(res, res, ct_buff, oper)
end

export Operator, geteval_at, ntt, ntt_to!, intt, intt_to!, divide_by_P, divide_by_P_to!, scale, scale_to!,
    divide_and_round, divide_and_round_to!, decompose, hoisted_gadgetprod_to!, gadgetprod, gadgetprod_to!,
    relinearise, relinearise_to!, extprod, extprod_to!, hoisted_extprod_to!