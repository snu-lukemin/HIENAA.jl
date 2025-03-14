struct Operator
    param::RingParam
    evalP::Union{Missing,PolyEvaluatorRNS}
    evalQ::PolyEvaluatorRNS
    auxeval::PolyEvaluatorArb
    decer::Decomposer
    ct_buff::Vector{RLWE}
    tensor_buff::Tensor

    function Operator(param::RLWEParameters)::Operator
        ring_param, P, Q, dlen = param.ring_param, param.P, param.Q, param.dlen
        N = ring_param.N

        if ismissing(P)
            Qlen = length(Q)
            Qmoduli = Modulus.(Q)
            evalP, evalQ = missing, PolyEvaluatorRNS(ring_param, Qmoduli)
            auxeval = PolyEvaluatorArb(ring_param, Modulus(1 << 62))
            decer = Decomposer(Qmoduli, dlen)
            ct_buff = [RLWE(N, Qlen, isPQ=false) for _ = 1:2]
            tensor_buff = Tensor(N, Qlen, decer.glen + 2, isPQ=true)
        else
            Plen, Qlen = length(P), length(Q)
            Pmoduli, Qmoduli = Modulus.(P), Modulus.(Q)
            evalP, evalQ = PolyEvaluatorRNS(ring_param, Pmoduli), PolyEvaluatorRNS(ring_param, Qmoduli)
            auxeval = PolyEvaluatorArb(ring_param, Modulus(1 << 62))
            decer = Decomposer(Pmoduli, Qmoduli, dlen)
            ct_buff = [RLWE(N, Plen + Qlen, isPQ=true) for _ = 1:2]
            tensor_buff = Tensor(N, Plen + Qlen, decer.glen + 2, isPQ=true)
        end

        new(ring_param, evalP, evalQ, auxeval, decer, ct_buff, tensor_buff)
    end
end

function geteval_at(len::Int64, oper::Operator; isPQ::Bool=false, auxQ::UInt64=UInt64(0))::PolyEvaluatorRNS
    evalP, evalQ = oper.evalP, oper.evalQ
    if isPQ
        if ismissing(evalP)
            throw(DomainError("Special modulus is not defined for the operator."))
        end
        Plen, Qlen = length(evalP), length(evalQ)
        if Qlen + Plen < len || Plen > len
            throw(BoundsError(len, "The length of the input polynomial is too large."))
        end
        eval = auxQ == 0 ? vcat(evalP, evalQ[1:len-Plen]) : vcat(evalP, evalQ[1:len-Plen-1], PolyEvaluatorArb(oper.auxeval, Modulus(auxQ)))
    else
        Qlen = length(evalQ)
        if len > Qlen
            throw(BoundsError(len, "The length of the input polynomial is too large."))
        end
        eval = auxQ == 0 ? evalQ[1:len] : vcat(evalQ[1:len-1], PolyEvaluatorArb(oper.auxeval, Modulus(auxQ)))
    end

    eval
end

#============================================================================================#
############################## Plaintext Operations ##########################################
#============================================================================================#

ntt(x::PlainPoly, oper::Operator)::PlainPoly = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::PlainPoly, x::PlainPoly, oper::Operator)::Nothing = begin
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output plaintexts have different length."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper; isPQ=isPQ, auxQ=auxQ)
    ntt_to!(res.val, x.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
    res.scale[] = x.scale[]

    return nothing
end

intt(x::PlainPoly, oper::Operator)::PlainPoly = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::PlainPoly, x::PlainPoly, oper::Operator)::Nothing = begin
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output plaintexts have different length."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper; isPQ=isPQ, auxQ=auxQ)
    intt_to!(res.val, x.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ
    res.scale[] = x.scale[]

    return nothing
end

neg(x::PlainText, oper::Operator)::PlainText = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

neg_to!(res::T, x::T, oper::Operator) where {T<:PlainText} = begin
    if length(res) ≠ length(x)
        throw(DimensionMismatch("The plaintexts should have the same length."))
    end

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    neg_to!(res.val, x.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
    res.scale[] = x.scale[]

    return nothing
end::Nothing

add(x::PlainText, y::PlainText, oper::Operator)::PlainText = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator)::Nothing = begin
    if length(res) ≠ length(x) ≠ length(y)
        throw(DimensionMismatch("The plaintexts should have the same length."))
    end
    if x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The plaintexts should be all in PQ or not in PQ."))
    end
    if x.auxQ[] ≠ y.auxQ[]
        throw(DomainError("The plaintexts should have the same auxiliary modulus."))
    end

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    add_to!(res.val, x.val, y.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
    res.scale[] = (x.scale[] + y.scale[]) / 2

    return nothing
end

sub(x::PlainText, y::PlainText, oper::Operator)::PlainText = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator)::Nothing = begin
    if length(res) ≠ length(x) ≠ length(y)
        throw(DimensionMismatch("The plaintexts should have the same length."))
    end
    if x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The plaintexts should be all in PQ or not in PQ."))
    end
    if x.auxQ[] ≠ y.auxQ[]
        throw(DomainError("The plaintexts should have the same auxiliary modulus."))
    end

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    sub_to!(res.val, x.val, y.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
    res.scale[] = (x.scale[] + y.scale[]) / 2

    return nothing
end

mul(x::PlainText, y::PlainText, oper::Operator)::PlainText = begin
    res = similar(x)
    mul_to!(res, x, y, oper)
    res
end

mul_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The plaintexts should have the same length."))
    end
    if x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The plaintexts should be all in PQ or not in PQ."))
    end
    if x.auxQ[] ≠ y.auxQ[]
        throw(DomainError("The plaintexts should have the same auxiliary modulus."))
    end

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    mul_to!(res.val, x.val, y.val, eval)
    res.isPQ[] = x.isPQ[]
    res.auxQ[] = x.auxQ[]
    res.scale[] = x.scale[] * y.scale[]

    return nothing
end

muladd_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The plaintexts should have the same length."))
    end
    if !(res.isPQ[] == x.isPQ[] == y.isPQ[])
        throw(DomainError("The plaintexts should be all in PQ or not in PQ."))
    end
    if !(res.auxQ[] == x.auxQ[] == y.auxQ[])
        throw(DomainError("The plaintexts should have the same auxiliary modulus."))
    end

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    muladd_to!(res.val, x.val, y.val, eval)
    res.scale[] = (res.scale[] + x.scale[] * y.scale[]) / 2

    return nothing
end

mulsub_to!(res::PlainText, x::PlainText, y::PlainText, oper::Operator)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The plaintexts should have the same length."))
    end
    if !(res.isPQ[] == x.isPQ[] == y.isPQ[])
        throw(DomainError("The plaintexts should be all in PQ or not in PQ."))
    end
    if !(res.auxQ[] == x.auxQ[] == y.auxQ[])
        throw(DomainError("The plaintexts should have the same auxiliary modulus."))
    end

    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    mulsub_to!(res.val, x.val, y.val, eval)
    res.scale[] = (res.scale[] + x.scale[] * y.scale[]) / 2

    return nothing
end

to_big(x::PlainText, oper::Operator) = begin
    eval = geteval_at(length(x), oper; isPQ=x.isPQ[], auxQ=x.auxQ[])
    to_big(x.val, eval)
end

change_modulus(x::PlainText, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))::PlainText = begin
    res = similar(x)
    change_modulus_to!(res, x, len, oper, auxQ=auxQ)
    res
end

function change_modulus_to!(res::PlainConst, x::PlainConst, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))::Nothing
    isPQ = x.isPQ[]
    now_Qlen, now_auxQ = length(x), x.auxQ[]
    tar_Qlen, tar_auxQ = len, auxQ

    if isPQ
        if ismissing(oper.evalP)
            throw(DomainError("The parameter does not support PQ."))
        end
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
        @views basis_extend_to!(res.val.vals[1:tar_Qlen], x.val.vals[1:now_Qlen], be)
    else
        @views basis_extend_to!(res.val.vals[1:tar_Qlen], x.val.vals[1:now_Qlen], be)
        resize!(res, tar_Qlen)
    end

    # Set parameters.
    res.auxQ[] = tar_auxQ
    res.isPQ[] = x.isPQ[]
    res.scale[] = x.scale[]

    return nothing
end

function change_modulus_to!(res::PlainPoly, x::PlainPoly, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))::Nothing
    isPQ = x.isPQ[]
    now_Qlen, now_auxQ = length(x), x.auxQ[]
    tar_Qlen, tar_auxQ = len, auxQ

    if isPQ
        if ismissing(oper.evalP)
            throw(DomainError("The parameter does not support PQ."))
        end
        tar_Qlen += length(oper.evalP)
    end

    # Easy cases.
    if tar_Qlen < now_Qlen && tar_auxQ == 0
        resize!(res, tar_Qlen)
        for i = 1:tar_Qlen
            copy!(res.val.coeffs[i], x.val.coeffs[i])
        end

        res.auxQ[] = tar_auxQ
        res.isPQ[] = x.isPQ[]
        res.val.isntt[] = x.val.isntt[]

        return
    elseif now_Qlen ≥ tar_Qlen && now_auxQ == tar_auxQ
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
    basis_extend_to!(res.val.coeffs, buff.coeffs, be)

    # NTT, with optimisation.
    if buff2.isntt[]
        reuselen = min(now_Qlen, tar_Qlen) - 1
        for i = 1:reuselen
            copy!(res.val.coeffs[i], buff2.coeffs[i])
        end
        for i = reuselen+1:tar_Qlen
            ntt!(res.val.coeffs[i], eval_tarQ[i])
        end
    end
    res.val.isntt[] = buff2.isntt[]

    # Set parameters.
    res.auxQ[] = tar_auxQ
    res.isPQ[] = x.isPQ[]
    res.scale[] = x.scale[]

    return nothing
end

scale(x::PlainPoly, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0))::PlainPoly = begin
    res = similar(x)
    scale_to!(res, x, len, oper, auxQ=auxQ)
    res
end

function scale_to!(res::PlainPoly, x::PlainPoly, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0), isntt::Bool=true)::Nothing
    if x.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end

    newlen, newauxQ = len, auxQ
    oldlen, oldauxQ = length(x), x.auxQ[]
    if newlen > len
        throw(DomainError("The new length should be less than or equal to the original length."))
    end

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
            Rinv.vals[j] = Bmul(Rinv.vals[j], invmod(R[i].Q, T[j].Q), T[j])
            TinvS.vals[i] = Bmul(TinvS.vals[i], invmod(T[j].Q, R[i].Q), R[i])
        end

        for i = eachindex(R), j = Rlen+1:newlen
            TinvS.vals[i] = Bmul(TinvS.vals[i], RS[j].Q, R[i])
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
        simple_scale_to!(res.val.coeffs, buffT.coeffs, ss)
        res.val.isntt[] = false

        if isntt
            ntt!(res.val, evalRS)
            for i = 1:length(evalRS)-1
                !buffRS.isntt[] && ntt!(buffRS[i], evalRS[i])
                add_to!(res.val.coeffs[i], res.val.coeffs[i], buffRS[i], evalRS[i])
            end
        else
            for i = 1:length(evalRS)-1
                buffRS.isntt[] && intt!(buffRS[i], evalRS[i])
                add_to!(res.val.coeffs[i], res.val.coeffs[i], buffRS[i], evalRS[i])
            end
        end

        res.scale[] = x.scale[] * (prod(evalRS.moduli) // prod(evalRT.moduli))
    else
        # Define Moduli.
        # R -> T
        evalR = geteval_at(oldlen, oper, auxQ=oldauxQ)
        evalT = geteval_at(newlen, oper, auxQ=newauxQ)

        ss = SimpleScaler(evalR.moduli, evalT.moduli)
        buff = oper.tensor_buff[1, 1:oldlen]

        copy!(buff, x.val)
        buff.isntt[] && intt!(buff, evalR)
        simple_scale_to!((@view res.val.coeffs[1:newlen]), buff.coeffs, ss)
        resize!(res, newlen)
        res.val.isntt[] = false
        isntt && ntt!(res.val, evalT)

        res.scale[] = x.scale[] * (prod(evalR.moduli) // prod(evalT.moduli))
    end

    res.isPQ[] = false
    res.auxQ[] = newauxQ

    return nothing
end

#============================================================================================#
############################## RLWE Operations ###############################################
#============================================================================================#

ntt(x::RLWE, oper::Operator)::RLWE = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

function ntt_to!(res::RLWE, x::RLWE, oper::Operator)::Nothing
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    ntt_to!(res.b, x.b, eval)
    ntt_to!(res.a, x.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

intt(x::RLWE, oper::Operator)::RLWE = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

function intt_to!(res::RLWE, x::RLWE, oper::Operator)::Nothing
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    intt_to!(res.b, x.b, eval)
    intt_to!(res.a, x.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

neg(x::RLWE, oper::Operator)::RLWE = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::RLWE, x::RLWE, oper::Operator)::Nothing
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    neg_to!(res.b, x.b, eval)
    neg_to!(res.a, x.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

add(x::T, y::S, oper::Operator) where {T,S<:Union{RLWE,PlainText}} = begin
    res = T == RLWE ? similar(x) : similar(y)
    add_to!(res, x, y, oper)
    res
end

function add_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    add_to!(res.b, x.val, y.b, eval)
    copy!(res.a, y.a)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

add_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator)::Nothing = begin
    add_to!(res, y, x, oper)
    return nothing
end

function add_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    add_to!(res.b, x.b, y.b, eval)
    add_to!(res.a, x.a, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

sub(x::T, y::S, oper::Operator) where {T,S<:Union{RLWE,PlainText}} = begin
    res = T == RLWE ? similar(x) : similar(y)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator)::Nothing
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    sub_to!(res.b, x.b, y.val, eval)
    copy!(res.a, x.a)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

function sub_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)::Nothing
    neg_to!(res, y, oper)
    add_to!(res, x, res, oper)
    return nothing
end

function sub_to!(res::RLWE, x::RLWE, y::RLWE, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    sub_to!(res.b, x.b, y.b, eval)
    sub_to!(res.a, x.a, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

mul(x::T, y::S, oper::Operator) where {T,S<:Union{RLWE,PlainText}} = begin
    res = T == RLWE ? similar(x) : similar(y)
    mul_to!(res, x, y, oper)
    res
end

function mul_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)::Nothing
    if !(length(x) == length(res) == length(y))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    mul_to!(res.b, x.val, y.b, eval)
    mul_to!(res.a, x.val, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

mul_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator)::Nothing = begin
    mul_to!(res, y, x, oper)
    return nothing
end

function muladd_to!(res::RLWE, x::PlainText, y::PlainText, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if !(x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[])
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    muladd_to!(res.b, x.val, y.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

function muladd_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if !(x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[])
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    muladd_to!(res.b, x.val, y.b, eval)
    muladd_to!(res.a, x.val, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

muladd_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator)::Nothing = begin
    muladd_to!(res, y, x, oper)
    return nothing
end

function mulsub_to!(res::RLWE, x::PlainText, y::PlainText, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if !(x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[])
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    mulsub_to!(res.b, x.val, y.val, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

function mulsub_to!(res::RLWE, x::PlainText, y::RLWE, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different length."))
    end
    if !(x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[])
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    mulsub_to!(res.b, x.val, y.b, eval)
    mulsub_to!(res.a, x.val, y.a, eval)

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

mulsub_to!(res::RLWE, x::RLWE, y::PlainText, oper::Operator)::Nothing = begin
    mulsub_to!(res, y, x, oper)
    return nothing
end

divide_by_P(x::RLWE, oper::Operator; isntt::Bool=true)::RLWE = begin
    res = similar(x)
    divide_by_P_to!(res, x, oper, isntt=isntt)
    res
end

function divide_by_P_to!(res::RLWE, x::RLWE, oper::Operator; isntt::Bool=true)::Nothing
    if !x.isPQ[]
        throw(DomainError("The input ciphertext should be in PQ."))
    end

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
        Pinv.vals[j] = Bmul(Pinv.vals[j], invmod(P[i].Q, Q[j].Q), Q[j])
        Qinv.vals[i] = Bmul(Qinv.vals[i], invmod(Q[j].Q, P[i].Q), P[i])
    end

    # scale b.
    copy!(buff, x.b)
    buffP.isntt[] = x.b.isntt[]
    buffQ.isntt[] = x.b.isntt[]

    mul_to!(buffP, Qinv, buffP, evalP)
    mul_to!(buffQ, Pinv, buffQ, evalQ)
    buffP.isntt[] && intt!(buffP, evalP)

    resize!(res.b, len)
    simple_scale_to!(res.b.coeffs, buffP.coeffs, ss)
    res.b.isntt[] = false

    if isntt
        ntt!(res.b, evalQ)
        !buffQ.isntt[] && ntt!(buffQ, evalQ)
    else
        buffQ.isntt[] && intt!(buffQ, evalQ)
    end
    add_to!(res.b, res.b, buffQ, evalQ)

    # scale a.
    copy!(buff, x.a)
    buffP.isntt[] = x.a.isntt[]
    buffQ.isntt[] = x.a.isntt[]

    mul_to!(buffP, Qinv, buffP, evalP)
    mul_to!(buffQ, Pinv, buffQ, evalQ)
    buffP.isntt[] && intt!(buffP, evalP)

    resize!(res.a, len)
    simple_scale_to!(res.a.coeffs, buffP.coeffs, ss)
    res.a.isntt[] = false

    if isntt
        ntt!(res.a, evalQ)
        !buffQ.isntt[] && ntt!(buffQ, evalQ)
    else
        buffQ.isntt[] && intt!(buffQ, evalQ)
    end
    add_to!(res.a, res.a, buffQ, evalQ)

    res.isPQ[] = false
    res.auxQ[] = auxQ

    return nothing
end

scale(x::RLWE, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0), isntt::Bool=true)::RLWE = begin
    res = similar(x)
    scale_to!(res, x, len, oper, auxQ=auxQ, isntt=isntt)
    res
end

function scale_to!(res::RLWE, x::RLWE, len::Int64, oper::Operator; auxQ::UInt64=UInt64(0), isntt::Bool=true)::Nothing
    if x.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end

    newlen, newauxQ = len, auxQ
    oldlen, oldauxQ = length(x), x.auxQ[]
    if newlen > len
        throw(DomainError("The new length should be less than or equal to the original length."))
    end

    if oldlen == newlen && oldauxQ == newauxQ
        copy!(res, x)
        
        return nothing
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
            Rinv.vals[j] = Bmul(Rinv.vals[j], invmod(R[i].Q, T[j].Q), T[j])
            TinvS.vals[i] = Bmul(TinvS.vals[i], invmod(T[j].Q, R[i].Q), R[i])
        end

        for i = eachindex(R), j = Rlen+1:newlen
            TinvS.vals[i] = Bmul(TinvS.vals[i], RS[j].Q, R[i])
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
        simple_scale_to!(res.b.coeffs, buffT.coeffs, ss)
        res.b.isntt[] = false

        if isntt
            ntt!(res.b, evalRS)
            for i = 1:length(evalRS)-1
                !buffRS.isntt[] && ntt!(buffRS[i], evalRS[i])
                add_to!(res.b.coeffs[i], res.b.coeffs[i], buffRS[i], evalRS[i])
            end
        else
            for i = 1:length(evalRS)-1
                buffRS.isntt[] && intt!(buffRS[i], evalRS[i])
                add_to!(res.b.coeffs[i], res.b.coeffs[i], buffRS[i], evalRS[i])
            end
        end

        # scale a.
        copy!(buffRT, x.a)
        buffR.isntt[] = x.a.isntt[]
        buffT.isntt[] = x.a.isntt[]
        buffRS.isntt[] = x.a.isntt[]
        resize!(res.a, newlen)

        mul_to!(buffR, TinvS, buffR, evalR)
        mul_to!(buffT, Rinv, buffT, evalT)
        buffT.isntt[] && intt!(buffT, evalT)
        simple_scale_to!(res.a.coeffs, buffT.coeffs, ss)
        res.a.isntt[] = false

        if isntt
            ntt!(res.a, evalRS)
            for i = 1:length(evalRS)-1
                !buffRS.isntt[] && ntt!(buffRS[i], evalRS[i])
                add_to!(res.a.coeffs[i], res.a.coeffs[i], buffRS[i], evalRS[i])
            end
        else
            for i = 1:length(evalRS)-1
                buffRS.isntt[] && intt!(buffRS[i], evalRS[i])
                add_to!(res.a.coeffs[i], res.a.coeffs[i], buffRS[i], evalRS[i])
            end
        end
    else
        # Define Moduli.
        # R -> T
        evalR = geteval_at(oldlen, oper, auxQ=oldauxQ)
        evalT = geteval_at(newlen, oper, auxQ=newauxQ)

        ss = SimpleScaler(evalR.moduli, evalT.moduli)
        buff = oper.tensor_buff[1, 1:oldlen]

        copy!(buff, x.b)
        buff.isntt[] && intt!(buff, evalR)
        simple_scale_to!((@view res.b.coeffs[1:newlen]), buff.coeffs, ss)
        resize!(res.b, newlen)
        res.b.isntt[] = false
        isntt && ntt!(res.b, evalT)

        copy!(buff, x.a)
        buff.isntt[] && intt!(buff, evalR)
        simple_scale_to!((@view res.a.coeffs[1:newlen]), buff.coeffs, ss)
        resize!(res.a, newlen)
        res.a.isntt[] = false
        isntt && ntt!(res.a, evalT)
    end

    res.isPQ[] = false
    res.auxQ[] = newauxQ

    return nothing
end

#============================================================================================#
############################## Tensor Operations #############################################
#============================================================================================#

ntt(x::Tensor, oper::Operator)::Tensor = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::Tensor, x::Tensor, oper::Operator)::Nothing = begin
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        ntt_to!(res.val[i], x.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

intt(x::Tensor, oper::Operator)::Tensor = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::Tensor, x::Tensor, oper::Operator)::Nothing = begin
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        intt_to!(res.val[i], x.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

neg(x::Tensor, oper::Operator)::Tensor = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

function neg_to!(res::Tensor, x::Tensor, oper::Operator)::Nothing
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        neg_to!(res.val[i], x.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

add(x::Tensor, y::Tensor, oper::Operator)::Tensor = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

function add_to!(res::Tensor, x::Tensor, y::Tensor, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        add_to!(res.val[i], x.val[i], y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

sub(x::Tensor, y::Tensor, oper::Operator)::Tensor = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

function sub_to!(res::Tensor, x::Tensor, y::Tensor, oper::Operator)::Nothing
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        sub_to!(res.val[i], x.val[i], y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

mul(x::Tensor, y::PlainText, oper::Operator)::Tensor = begin
    res = similar(x)
    mul_to!(res, x, y, oper)
    res
end

mul(x::PlainText, y::Tensor, oper::Operator)::Tensor = mul(y, x, oper)

function mul_to!(res::Tensor, x::PlainText, y::Tensor, oper::Operator)::Nothing
    if !(length(x) == length(res) == length(y))
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end
    if x.auxQ[] ≠ y.auxQ[] || x.isPQ[] ≠ y.isPQ[]
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        mul_to!(res.val[i], x.val, y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

mul_to!(res::Tensor, x::Tensor, y::PlainText, oper::Operator)::Nothing = begin
    mul_to!(res, y, x, oper)
    return nothing
end

muladd_to!(res::Tensor, x::PlainText, y::Tensor, oper::Operator)::Nothing = begin
    if !(length(x) == length(res) == length(y))
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end
    if !(x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[])
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        muladd_to!(res.val[i], x.val, y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

muladd_to!(res::Tensor, x::Tensor, y::PlainText, oper::Operator)::Nothing = begin
    muladd_to!(res, y, x, oper)
    return nothing
end

mulsub_to!(res::Tensor, x::PlainText, y::Tensor, oper::Operator)::Nothing = begin
    if !(length(x) == length(y) == length(res))
        throw(DimensionMismatch("The input and output ciphertexts have different sizes."))
    end
    if !(x.auxQ[] == y.auxQ[] == res.auxQ[] && x.isPQ[] == y.isPQ[] == res.isPQ[])
        throw(DomainError("The input and output ciphertexts have different moduli."))
    end

    isPQ, auxQ = x.isPQ[], x.auxQ[]
    eval = geteval_at(length(x), oper, isPQ=isPQ, auxQ=auxQ)

    for i = 1:N
        mulsub_to!(res.val[i], x.val, y.val[i], eval)
    end

    res.isPQ[] = isPQ
    res.auxQ[] = auxQ

    return nothing
end

mulsub_to!(res::Tensor, x::Tensor, y::PlainText, oper::Operator)::Nothing = begin
    mulsub_to!(res, y, x, oper)
    return nothing
end

#================================================================================================#
##################################### RLEV Operations ############################################
#================================================================================================#

ntt(x::RLEV, oper::Operator)::RLEV = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::RLEV, x::RLEV, oper::Operator)::Nothing = begin
    if res.glen ≠ x.glen
        throw(DimensionMismatch("The length of input and output ciphertext should match."))
    end

    for i = 1:res.len
        ntt_to!(res.stack[i], x.stack[i], oper)
    end

    return nothing
end

intt(x::RLEV, oper::Operator)::RLEV = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::RLEV, x::RLEV, oper::Operator)::Nothing = begin
    if res.glen ≠ x.glen
        throw(DimensionMismatch("The length of input and output ciphertext should match."))
    end

    for i = 1:res.len
        intt_to!(res.stack[i], x.stack[i], oper)
    end

    return nothing
end

neg(x::RLEV, oper::Operator)::RLEV = begin
    res = similar(x)
    neg_to!(res, x, oper)
    res
end

neg_to!(res::RLEV, x::RLEV, oper::Operator)::Nothing = begin
    if res.glen ≠ x.glen
        throw(DimensionMismatch("The length of input and output ciphertext should match."))
    end

    for i = 1:res.len
        neg_to!(res.stack[i], x.stack[i], oper)
    end

    return nothing
end

add(x::RLEV, y::RLEV, oper::Operator)::RLEV = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add_to!(res::RLEV, x::RLEV, y::RLEV, oper::Operator)::Nothing = begin
    if !(res.glen == x.glen == y.glen)
        throw(DimensionMismatch("The length of input and output ciphertext should match."))
    end

    for i = 1:res.len
        add_to!(res.stack[i], x.stack[i], y.stack[i], oper.eval)
    end

    return nothing
end

sub(x::RLEV, y::RLEV, oper::Operator)::RLEV = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::RLEV, x::RLEV, y::RLEV, oper::Operator)::Nothing = begin
    if !(res.glen == x.glen == y.glen)
        throw(DimensionMismatch("The length of input and output ciphertext should match."))
    end

    for i = 1:res.len
        sub_to!(res.stack[i], x.stack[i], y.stack[i], oper.eval)
    end

    return nothing
end

#================================================================================================#
##################################### RGSW Operations ############################################
#================================================================================================#

ntt(x::RGSW, oper::Operator)::RGSW = begin
    res = similar(x)
    ntt_to!(res, x, oper)
    res
end

ntt_to!(res::RGSW, x::RGSW, oper::Operator)::Nothing = begin
    ntt_to!(res.basketb, x.basketb, oper)
    ntt_to!(res.basketa, x.basketa, oper)

    return nothing
end

intt(x::RGSW, oper::Operator)::RGSW = begin
    res = similar(x)
    intt_to!(res, x, oper)
    res
end

intt_to!(res::RGSW, x::RGSW, oper::Operator)::Nothing = begin
    intt_to!(res.basketb, x.basketb, oper)
    intt_to!(res.basketa, x.basketa, oper)

    return nothing
end

add(x::RGSW, y::RGSW, oper::Operator)::RGSW = begin
    res = similar(x)
    add_to!(res, x, y, oper)
    res
end

add_to!(res::RGSW, x::RGSW, y::RGSW, oper::Operator)::Nothing = begin
    add_to!(res.basketb, x.basketb, y.basketb, oper.eval)
    add_to!(res.basketa, x.basketa, y.basketa, oper.eval)

    return nothing
end

sub(x::RGSW, y::RGSW, oper::Operator)::RGSW = begin
    res = similar(x)
    sub_to!(res, x, y, oper)
    res
end

sub_to!(res::RGSW, x::RGSW, y::RGSW, oper::Operator)::Nothing = begin
    sub_to!(res.basketb, x.basketb, y.basketb, oper.eval)
    sub_to!(res.basketa, x.basketa, y.basketa, oper.eval)

    return nothing
end

#================================================================================================#
##################################### Decompositions #############################################
#================================================================================================#

decompose(x::PlainPoly, oper::Operator)::Tensor = begin
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

function decompose_to!(decx::Tensor, x::PlainPoly, oper::Operator)::Nothing
    if x.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end

    len, auxQ, decer = length(x), x.auxQ[], oper.decer
    buff = PlainPoly(oper.tensor_buff[end, 1:len])

    copy!(buff, x)
    buff.val.isntt[] && intt_to!(buff, buff, oper)
    auxQ ≠ 0 && scale_to!(buff, buff, len, oper, isntt=false)

    buff.isPQ[] = false
    buff.auxQ[] = auxQ

    decompose_to!(decx, buff, decer)

    return nothing
end

function hoisted_gadgetprod_to!(res::RLWE, decx::Tensor, ct::RLEV, oper::Operator; islazy::Bool=false)::Nothing
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
        if ismissing(oper.evalP)
            throw(DomainError("The ciphertext is not defined over PQ."))
        end

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
        scale_to!(res, res, len, oper, auxQ=auxQ)
        evalQ = geteval_at(len, oper, auxQ=auxQ)
    end

    !res.b.isntt[] && ntt!(res.b, evalQ)
    !res.a.isntt[] && ntt!(res.a, evalQ)

    return nothing
end

gadgetprod(x::PlainPoly, ct::RLEV, oper::Operator)::RLWE = begin
    res = similar(ct.stack[1])
    gadgetprod_to!(res, x, ct, oper)
    res
end

function gadgetprod_to!(res::RLWE, x::PlainPoly, ct::RLEV, oper::Operator; islazy::Bool=false)::Nothing
    len, decer = length(x), oper.decer

    if ismissing(oper.evalP)
        dlen = ceil(Int64, len / decer.dlen)
        decbuff = oper.tensor_buff[1:dlen, 1:len]
    else
        dlen, Plen = ceil(Int64, len / decer.dlen), length(oper.evalP)
        decbuff = oper.tensor_buff[1:dlen, 1:len+Plen]
    end

    decompose_to!(decbuff, x, oper)
    hoisted_gadgetprod_to!(res, decbuff, ct, oper, islazy=islazy)

    return nothing
end

relinearise(ct::Tensor, rlk::RLEV, oper::Operator)::RLWE = begin
    res = RLWE(similar(ct.vals[1]), similar(ct.vals[1]))
    relinearise_to!(res, ct, rlk, oper)
    res
end

function relinearise_to!(res::RLWE, ct::Tensor, rlk::RLEV, oper::Operator)::Nothing
    if ct.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end
    if length(ct.vals) ≠ 3
        throw(DimensionMismatch("The input tensor length should be 3."))
    end

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

    return nothing
end

keyswitch(ct::RLWE, ksk::RLEV, oper::Operator)::RLWE = begin
    res = similar(ct)
    keyswitch_to!(res, ct, ksk, oper)
    res
end

function keyswitch_to!(res::RLWE, ct::RLWE, ksk::RLEV, oper::Operator)::Nothing
    if ct.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = geteval_at(len, oper, auxQ=auxQ)
    buff = oper.tensor_buff[end-1, 1:len]

    copy!(buff, ct.b)
    !buff.isntt[] && ntt!(buff, evalQ)

    gadgetprod_to!(res, PlainPoly(ct.a, auxQ=auxQ), ksk, oper)

    add_to!(res.b, res.b, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]

    return nothing
end

hoisted_keyswitch(adec::Tensor, ct::RLWE, ksk::RLEV, oper::Operator)::RLWE = begin
    res = similar(ct)
    hoisted_keyswitch_to!(res, adec, ct, ksk, oper)
    res
end

function hoisted_keyswitch_to!(res::RLWE, adec::Tensor, ct::RLWE, ksk::RLEV, oper::Operator)::Nothing
    if ct.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end

    len, auxQ = length(ct), ct.auxQ[]
    evalQ = geteval_at(len, oper, auxQ=auxQ)
    buff = oper.tensor_buff[end-1, 1:len]

    copy!(buff, ct.b)
    !buff.isntt[] && ntt!(buff, evalQ)

    hoisted_gadgetprod_to!(res, adec, ksk, oper)

    add_to!(res.b, res.b, buff, evalQ)

    res.isPQ[] = false
    res.auxQ[] = ct.auxQ[]

    return nothing
end

automorphism(x::RLWE, idx::Integer, atk::RLEV, oper::Operator)::RLWE = begin
    res = deepcopy(x)
    automorphism_to!(res, x, idx, atk, oper)
    res
end

automorphism!(x::RLWE, idx::Integer, atk::RLEV, oper::Operator)::Nothing = begin
    automorphism_to!(x, x, idx, atk, oper)
    return nothing
end

function automorphism_to!(res::RLWE, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)::Nothing
    eval = geteval_at(length(ct), oper, isPQ=ct.isPQ[], auxQ=ct.auxQ[])

    keyswitch_to!(res, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)

    return nothing
end

hoisted_automorphism(adec::Tensor, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)::RLWE = begin
    res = similar(ct)
    hoisted_automorphism_to!(res, adec, ct, idx, atk, oper)
    res
end

function hoisted_automorphism_to!(res::RLWE, adec::Tensor, ct::RLWE, idx::Integer, atk::RLEV, oper::Operator)::Nothing
    eval = geteval_at(length(ct), oper, isPQ=ct.isPQ[], auxQ=ct.auxQ[])

    hoisted_keyswitch_to!(res, adec, ct, atk, oper)
    automorphism!(res.b, idx, eval)
    automorphism!(res.a, idx, eval)

    return nothing
end

#======================================================================================#

function extprod(ct::RLWE, rgsw::RGSW, oper::Operator)::RLWE
    res = similar(ct)
    extprod_to!(res, ct, rgsw, oper)
    res
end

function extprod_to!(res::RLWE, ct::RLWE, rgsw::RGSW, oper::Operator)::Nothing
    if ct.isPQ[]
        throw(DomainError("The input ciphertext should not be in PQ."))
    end

    len, auxQ = length(ct), ct.auxQ[]
    buff = oper.tensor_buff[end-1, 1:len]
    ct_buff = oper.ct_buff[1][1:len]

    copy!(buff, ct.a)
    gadgetprod_to!(res, PlainPoly(ct.b, auxQ=auxQ), rgsw.basketb, oper)
    gadgetprod_to!(ct_buff, PlainPoly(buff, auxQ=auxQ), rgsw.basketa, oper)
    add_to!(res, res, ct_buff, oper)

    return nothing
end

@views function hoisted_extprod_to!(res::RLWE, ctdec::Tensor, rgsw::RGSW, oper::Operator)::Nothing
    len = length(ctdec)
    !ismissing(oper.evalP) && (len -= length(oper.evalP))
    ct_buff, glen = oper.ct_buff[1][1:len], length(ctdec) >> 1

    hoisted_gadgetprod_to!(res, ctdec[1:glen], rgsw.basketb, oper)
    hoisted_gadgetprod_to!(ct_buff, ctdec[glen+1:end], rgsw.basketa, oper)
    add_to!(res, res, ct_buff, oper)

    return nothing
end
