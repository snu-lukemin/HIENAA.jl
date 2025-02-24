"""
SKEncryptor is a struct for secret-key encryption.
"""
struct SKEncryptor
    key::RLWEkey
    keyPQ::RLWEkeyPQ
    usampler::UniformSampler
    gsampler::CDTSampler
    oper::Operator

    SKEncryptor(key::RLWEkey, σ::Real, oper::Operator) = begin
        evalP, evalQ = oper.evalP, oper.evalQ

        if ismissing(evalP)
            keyPQ = RLWEkeyPQ(key, evalQ)
            ntt!(keyPQ, evalQ)
        else
            evalPQ = vcat(evalP, evalQ)

            keyPQ = RLWEkeyPQ(key, evalPQ)
            ntt!(keyPQ, evalPQ)
        end

        new(key, keyPQ, UniformSampler(), CDTSampler(0.0, σ), oper)
    end
end

_getkey_at(len::Int64, entor::SKEncryptor; isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    evalP, evalQ = entor.oper.evalP, entor.oper.evalQ

    if isPQ
        @assert !ismissing(evalP) "The key is not defined over R_PQ."
        Plen, Qlen = length(evalP), length(evalQ)
        @assert 0 < len ≤ Plen + Qlen "The length of input RLWE ciphertext does not match the parameters."
        if auxQ == 0
            key = entor.keyPQ[1:len]
        else
            key = entor.keyPQ[1:len]
            key[len] = entor.oper.tensor_buff[1, 1]
            Q = Modulus(auxQ)
            @inbounds for i = 1:key.N
                key[len][i] = _Bred(entor.key.coeffs[i], Q)
            end
        end
    else
        Qlen = length(evalQ)
        @assert 0 < len ≤ Qlen "The length of input RLWE ciphertext does not match the parameters."
        if auxQ == 0
            if ismissing(evalP)
                key = entor.keyPQ[1:len]
            else
                Plen = length(evalP)
                key = entor.keyPQ[Plen+1:Plen+len]
            end
        else
            if ismissing(evalP)
                key = entor.keyPQ[1:len]
                key[len] = entor.oper.tensor_buff[1, 1]
                Q = Modulus(auxQ)
                @inbounds for i = 1:key.N
                    key[len][i] = _Bred(entor.key.coeffs[i], Q)
                end
            else
                Plen = length(evalP)
                key = entor.keyPQ[Plen+1:Plen+len]
                key[len] = entor.oper.tensor_buff[1, 1]
                Q = Modulus(auxQ)
                @inbounds for i = 1:key.N
                    key[len][i] = _Bred(entor.key.coeffs[i], Q)
                end
            end
        end
    end

    key
end

function rlwe_sample(entor::SKEncryptor, Qlen::Int64=typemax(Int64); isPQ::Bool=false, auxQ::UInt64=UInt64(0))
    oper, N = entor.oper, entor.keyPQ.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_sample_to!(res, entor)
    res
end

function rlwe_sample_to!(res::RLWE, entor::SKEncryptor)
    us, gs = entor.usampler, entor.gsampler
    N, len, isPQ, auxQ = res.a.N, length(res.a), res.isPQ[], res.auxQ[]

    key = _getkey_at(len, entor, isPQ=isPQ, auxQ=auxQ)
    eval = geteval_at(len, entor.oper, isPQ=isPQ, auxQ=auxQ)

    b, a = res.b, res.a
    b.isntt[] = false
    a.isntt[] = true

    uniform_random_to!(us, a, eval)

    @inbounds for j = 1:N
        ei = sample(gs)
        for i = 1:len
            b.coeffs[i][j] = _Bred(ei, eval[i])
        end
    end

    ntt!(b, eval)
    mulsub_to!(b, a, key, eval)
end

public_keygen(entor::SKEncryptor) = rlwe_sample(entor, isPQ=!ismissing(entor.oper.evalP))

"""
PKEncryptor is a struct for public key encryptions.
"""
struct PKEncryptor
    pk::RLWE
    gsampler::CDTSampler
    rgsampler::RGSampler
    oper::Operator

    PKEncryptor(pk::RLWE, σ::Real, τ::Real, oper::Operator) = begin
        evalP, evalQ = oper.evalP, oper.evalQ

        if ismissing(evalP)
            len = length(evalQ)

            @assert !pk.isPQ[] && length(pk) == len "The public key length does not match the parameters."
        else
            len = length(evalP) + length(evalQ)

            @assert pk.isPQ[] && length(pk) == len "The public key length does not match the parameters."
        end

        new(pk, CDTSampler(0.0, σ), RGSampler(τ), oper)
    end
end

function rlwe_sample(entor::PKEncryptor, Qlen::Int64=typemax(Int64); isPQ::Bool=false, auxQ::UInt64=UInt64(0))
    oper, N = entor.oper, entor.oper.param.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_sample_to!(res, entor)
    res
end

function rlwe_sample_to!(res::RLWE, entor::PKEncryptor)
    pk, gs, rgs, evalP = entor.pk, entor.gsampler, entor.rgsampler, entor.oper.evalP
    N, len, isPQ, auxQ = res.a.N, length(res.a), res.isPQ[], res.auxQ[]

    if isPQ
        pk = entor.pk[1:len]
    else
        if ismissing(evalP)
            pk = entor.pk[1:len]
        else
            Plen = length(evalP)
            pk = entor.pk[Plen+1:Plen+len]
        end
    end
    eval = geteval_at(len, entor.oper, isPQ=isPQ)

    buff = entor.oper.tensor_buff[1, 1:len]
    b, a = res.b, res.a

    buff.isntt[] = false
    @inbounds for j = 1:N
        rj = sample(gs)
        for i = 1:len
            buff.coeffs[i][j] = _Bred(rj, eval[i])
        end
    end

    ntt!(buff, eval)
    mul_to!(b, buff, pk.b, eval)
    mul_to!(a, buff, pk.a, eval)

    buff.isntt[] = false
    @inbounds for j = 1:N
        e1j = sample(gs)
        for i = 1:len
            buff.coeffs[i][j] = _Bred(e1j, eval[i])
        end
    end

    ntt!(buff, eval)
    add_to!(a, buff, a, eval)

    buff.isntt[] = false
    @inbounds for j = 1:N
        e0j = sample(rgs)
        for i = 1:len
            buff.coeffs[i][j] = _Bred(e0j, eval[i])
        end
    end

    ntt!(buff, eval)
    add_to!(b, buff, b, eval)

    if auxQ ≠ 0
        intt!(b, eval)
        intt!(a, eval)

        oldQ = eval.moduli
        newQ = vcat(oldQ[1:end-1], Modulus(auxQ))

        ss = SimpleScaler(oldQ, newQ)
        eval = geteval_at(len, entor.oper, isPQ=isPQ, auxQ=auxQ)

        simple_scale!(res.b.coeffs, res.b.coeffs, ss)
        res.b.isntt[] = false

        simple_scale!(res.a.coeffs, res.a.coeffs, ss)
        res.a.isntt[] = false
    end
end

"""
Encryptor is a struct for encryption and decryption of RLWE-based ciphertexts.
"""
Encryptor = Union{SKEncryptor,PKEncryptor}

(::Type{Encryptor})(key::RLWEkey, σ::Real, oper::Operator) = SKEncryptor(key, σ, oper)
(::Type{Encryptor})(pk::RLWE, σ::Real, τ::Real, oper::Operator) = PKEncryptor(pk, σ, τ, oper)

"""
rlwe_encrypt is a function to encrypt a plaintext into RLWE ciphertext.
"""
rlwe_encrypt(m::Union{PlainConst,PlainPoly}, entor::Encryptor) = begin
    oper, N = entor.oper, entor.oper.param.N
    m.isPQ[] && @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."

    res = RLWE(N, length(m))
    rlwe_encrypt_to!(res, m, entor)
    res
end

rlwe_encrypt_to!(res::RLWE, m::PlainConst, entor::Encryptor) = begin
    @assert length(m) == length(res) "The length of input plaintext does not match the parameters."

    res.isPQ[], res.auxQ[] = m.isPQ[], m.auxQ[]
    rlwe_sample_to!(res, entor)
    eval = geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = deepcopy(m.val)

    add_to!(res.b, res.b, buff, eval)
    res
end

rlwe_encrypt_to!(res::RLWE, m::PlainPoly, entor::Encryptor) = begin
    @assert length(m) == length(res) "The length of input plaintext does not match the parameters."

    res.isPQ[], res.auxQ[] = m.isPQ[], m.auxQ[]
    rlwe_sample_to!(res, entor)
    eval = geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = entor.oper.tensor_buff[1, 1:len]

    copy!(buff, m.val)
    !buff.isntt[] && ntt!(buff, eval)

    add_to!(res.b, res.b, buff, entor.eval)
    res
end

rlwe_encrypt_a(m::Union{PlainConst,PlainPoly}, entor::Encryptor) = begin
    oper, N = entor.oper, entor.oper.param.N
    if isPQ
        @assert !ismissing(oper.evalP) "The key is not defined over R_PQ."
        len = length(oper.evalP) + min(Qlen, length(oper.evalQ))
    else
        len = min(Qlen, length(oper.evalQ))
    end
    res = RLWE(N, len, isPQ=isPQ, auxQ=auxQ)
    rlwe_encrypt_a_to!(res, m, entor)
    res
end

rlwe_encrypt_a_to!(res::RLWE, m::PlainConst, entor::Encryptor) = begin
    @assert length(m) == length(res) "The length of input plaintext does not match the parameters."

    res.isPQ[], res.auxQ[] = m.isPQ[], m.auxQ[]
    rlwe_sample_to!(res, entor)
    eval = geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = deepcopy(m.val)

    add_to!(res.a, res.a, buff, eval)
    res
end

rlwe_encrypt_a_to!(res::RLWE, m::PlainPoly, entor::Encryptor) = begin
    @assert length(m) == length(res) "The length of input plaintext does not match the parameters."

    res.isPQ[], res.auxQ[] = m.isPQ[], m.auxQ[]
    rlwe_sample_to!(res, entor)
    eval = geteval_at(length(res), entor.oper, isPQ=res.isPQ[], auxQ=res.auxQ[])
    buff = entor.oper.tensor_buff[1, 1:len]

    copy!(buff, m.val)
    !buff.isntt[] && ntt!(buff, eval)

    add_to!(res.a, res.a, buff, entor.eval)
    res
end

function phase(ct::RLWE, entor::SKEncryptor)
    res = PlainPoly(similar(ct.a))
    phase_to!(res, ct, entor)
    res
end

function phase_to!(res::PlainPoly, ct::RLWE, entor::SKEncryptor)
    @assert length(res) == length(ct) "The length of the output plaintext and input ciphertext should match."
    len, isPQ, auxQ = length(ct), ct.isPQ[], ct.auxQ[]
    buff = entor.oper.tensor_buff[end, 1:len]

    key = _getkey_at(len, entor, isPQ=isPQ, auxQ=auxQ)
    eval = geteval_at(len, entor.oper, isPQ=isPQ, auxQ=auxQ)

    copy!(res.val, ct.b)
    copy!(buff, ct.a)

    !res.val.isntt[] && ntt!(res.val, eval)
    !buff.isntt[] && ntt!(buff, eval)

    muladd_to!(res.val, buff, key, eval)

    intt!(res.val, eval)
    res.isPQ[], res.auxQ[] = isPQ, auxQ
end

#=====================================================================================================#

rlev_encrypt(m::Union{PlainConst,PlainPoly}, entor::Encryptor) = begin
    N, evalP, decer = entor.oper.param.N, entor.oper.evalP, entor.oper.decer

    @assert !m.isPQ[] "The messages cannot be defined over PQ."

    mlen = length(m)
    glen = ceil(Int64, mlen / decer.dlen)

    if ismissing(evalP)
        res = RLEV(N, mlen, glen)
    else
        res = RLEV(N, mlen + length(evalP), glen)
    end

    rlev_encrypt_to!(res, m, entor)
    res
end

function rlev_encrypt_to!(res::RLEV, m::PlainConst, entor::Encryptor)
    decer = entor.oper.decer

    @assert !m.isPQ[] "The messages cannot be defined over PQ."
    @assert m.auxQ[] == 0 "RLEV encryptions cannot be encrypted with auxiliary modulus."

    mlen = length(m)
    eval = geteval_at(mlen, entor.oper)
    glen = ceil(Int64, mlen / decer.dlen)
    if !ismissing(entor.oper.evalP)
        Plen = length(entor.oper.evalP)

        @inbounds for i = 1:glen
            resize!(res.stack[i], mlen + Plen)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)
            for j = 1:mlen
                _add_to!(res.stack[i].b[Plen+j], res.stack[i].b[Plen+j], _mul(m.val[j], decer.gvec[i][Plen+j], eval[j]), true, eval[j])
            end
        end
    else
        @inbounds for i = 1:glen
            resize!(res.stack[i], mlen)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            for j = 1:mlen
                _add_to!(res.stack[i].b[j], res.stack[i].b[j], _mul(m.val[j], decer.gvec[i][j], eval[j]), true, eval[j])
            end
        end
    end
end

function rlev_encrypt_to!(res::RLEV, m::PlainPoly, entor::Encryptor)
    decer = entor.oper.decer

    @assert !m.isPQ[] "The messages cannot be defined over PQ."
    @assert m.auxQ[] == 0 "RLEV encryptions cannot be encrypted with auxiliary modulus."

    mlen = length(m)
    eval = geteval_at(mlen, entor.oper)
    glen = ceil(Int64, mlen / decer.dlen)
    if !ismissing(entor.oper.evalP)
        Plen = length(entor.oper.evalP)

        buff = entor.oper.tensor_buff[2, 1:mlen]
        copy!(buff, m.val)
        !buff.isntt[] && ntt!(buff, eval)

        for i = 1:glen
            resize!(res.stack[i], mlen + Plen)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            for j = 1:mlen
                _muladd_to!(res.stack[i].b[Plen+j], decer.gvec[i][Plen+j], buff[j], eval[j])
            end
        end
    else
        buff = entor.oper.tensor_buff[2, 1:mlen]
        copy!(buff, m.val)
        !buff.isntt[] && ntt!(buff, eval)

        @inbounds for i = 1:glen
            resize!(res.stack[i], mlen)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            for j = 1:mlen
                _muladd_to!(res.stack[i].b[j], decer.gvec[i][j], buff[j], eval[j])
            end
        end
    end
end

rlev_encrypt_a(m::Union{PlainConst,PlainPoly}, entor::Encryptor) = begin
    N, evalP, decer = entor.oper.param.N, entor.oper.evalP, entor.oper.decer

    @assert !m.isPQ[] "The messages cannot be defined over PQ."

    mlen = length(m)
    glen = ceil(Int64, mlen / decer.dlen)

    if ismissing(evalP)
        res = RLEV(N, mlen, glen)
    else
        res = RLEV(N, mlen + length(evalP), glen)
    end

    rlev_encrypt_a_to!(res, m, entor)
    res
end

function rlev_encrypt_a_to!(res::RLEV, m::PlainConst, entor::Encryptor)
    decer = entor.oper.decer

    @assert !m.isPQ[] "The messages cannot be defined over PQ."
    @assert m.auxQ[] == 0 "RLEV encryptions cannot be encrypted with auxiliary modulus."

    mlen = length(m)
    eval = geteval_at(mlen, entor.oper)
    glen = ceil(Int64, mlen / decer.dlen)
    if !ismissing(entor.oper.evalP)
        Plen = length(entor.oper.evalP)

        @inbounds for i = 1:glen
            resize!(res.stack[i], mlen + Plen)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)
            for j = 1:mlen
                _add_to!(res.stack[i].a[Plen+j], res.stack[i].a[Plen+j], _mul(m.val[j], decer.gvec[i][Plen+j], eval[j]), true, eval[j])
            end
        end
    else
        @inbounds for i = 1:glen
            resize!(res.stack[i], mlen)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            for j = 1:mlen
                _add_to!(res.stack[i].a[j], res.stack[i].a[j], _mul(m.val[j], decer.gvec[i][j], eval[j]), true, eval[j])
            end
        end
    end
end

function rlev_encrypt_a_to!(res::RLEV, m::PlainPoly, entor::Encryptor)
    decer = entor.oper.decer

    @assert !m.isPQ[] "The messages cannot be defined over PQ."
    @assert m.auxQ[] == 0 "RLEV encryptions cannot be encrypted with auxiliary modulus."

    mlen = length(m)
    eval = geteval_at(mlen, entor.oper)
    glen = ceil(Int64, mlen / decer.dlen)
    if !ismissing(entor.oper.evalP)
        Plen = length(entor.oper.evalP)

        buff = entor.oper.tensor_buff[2, 1:mlen]
        copy!(buff, m.val)
        !buff.isntt[] && ntt!(buff, eval)

        for i = 1:glen
            resize!(res.stack[i], mlen + Plen)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            for j = 1:mlen
                _muladd_to!(res.stack[i].a[Plen+j], decer.gvec[i][Plen+j], buff[j], eval[j])
            end
        end
    else
        buff = entor.oper.tensor_buff[2, 1:mlen]
        copy!(buff, m.val)
        !buff.isntt[] && ntt!(buff, eval)

        @inbounds for i = 1:glen
            resize!(res.stack[i], mlen)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            for j = 1:mlen
                _muladd_to!(res.stack[i].a[j], decer.gvec[i][j], buff[j], eval[j])
            end
        end
    end
end

#=====================================================================================================#

function rgsw_encrypt(m::Union{PlainConst,PlainPoly}, entor::Encryptor)
    basketb = rlev_encrypt(m, entor)
    basketa = rlev_encrypt_a(m, entor)
    RGSW(basketb, basketa)
end

function rgsw_encrypt_to!(res::RGSW, m::Union{PlainConst,PlainPoly}, entor::Encryptor)
    rlev_encrypt_to!(res.basketb, m, entor)
    rlev_encrypt_a_to!(res.basketa, m, entor)
end

#=====================================================================================================#

function relin_keygen(entor::SKEncryptor)
    evalQ, evalP, decer = entor.oper.evalQ, entor.oper.evalP, entor.oper.decer

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)

        res = RLEV(entor.oper.param.N, len, decer.glen)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], entor.keyPQ, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)

        res = RLEV(entor.oper.param.N, len, decer.glen)

        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)

            muladd_to!(res.stack[i].a, decer.gvec[i], entor.keyPQ, eval)
        end
    end

    res
end

function automorphism_keygen(idx::Int64, entor::SKEncryptor)
    buff, evalQ, evalP, decer, param = entor.oper.tensor_buff.vals[2], entor.oper.evalQ, entor.oper.evalP, entor.oper.decer, entor.oper.param

    if ismissing(evalP)
        eval, len = evalQ, length(evalQ)
        res = RLEV(param.N, len, decer.glen)

        copy!(buff, entor.keyPQ)
        invidx = invmod(idx, param.m)
        automorphism!(buff, idx, eval)
        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = false
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)
            muladd_to!(res.stack[i].b, decer.gvec[i], buff, eval)
            automorphism!(res.stack[i].b, invidx, eval)
            automorphism!(res.stack[i].a, invidx, eval)
        end
    else
        eval, len = vcat(evalP, evalQ), length(evalP) + length(evalQ)
        res = RLEV(param.N, len, decer.glen)

        copy!(buff, entor.keyPQ)
        invidx = invmod(idx, param.m)
        automorphism!(buff, idx, eval)
        @inbounds for i = 1:res.glen
            resize!(res.stack[i], len)
            res.stack[i].isPQ[] = true
            res.stack[i].auxQ[] = 0

            rlwe_sample_to!(res.stack[i], entor)
            muladd_to!(res.stack[i].b, decer.gvec[i], buff, eval)
            automorphism!(res.stack[i].b, invidx, eval)
            automorphism!(res.stack[i].a, invidx, eval)
        end
    end

    res
end

export SKEncryptor, PKEncryptor, Encryptor, rlwe_sample, rlwe_sample_to!, rlwe_encrypt, rlwe_encrypt_to!, rlwe_encrypt_a, rlwe_encrypt_a_to!, 
    phase, phase_to!, rlev_encrypt, rlev_encrypt_to!, rlev_encrypt_a, rlev_encrypt_a_to!, rgsw_encrypt, rgsw_encrypt_to!