"""
    BGVScheme(sketch::BGVParamSketch)
    BGVScheme(param::BGVParameters)
    BGVScheme(oper::BGVOperator)

BGVScheme is a struct for the BGV homomorphic encryption scheme. It contains the following fields: oper, entor, rlk, and atk.
- oper: BGVOperator
- entor: Union{Missing,Encryptor}
- rlk: Union{Missing,RLEV}
- atk: Dict{Int64,RLEV}
"""
mutable struct BGVScheme
    oper::BGVOperator
    entor::Union{Missing,Encryptor}
    rlk::Union{Missing,RLEV}
    atk::Dict{Int64,RLEV}

    BGVScheme(sketch::BGVParamSketch) = BGVScheme(BGVParameters(sketch))

    function BGVScheme(param::BGVParameters)
        oper = BGVOperator(param)

        new(oper, missing, missing, Dict{Int64,RLEV}())
    end

    BGVScheme(oper::BGVOperator) = new(oper, missing, missing, Dict{Int64,RLEV}())
end

"""
    encrypt(msg::PlainText, scheme::HEScheme)

Encrypts a plaintext message `msg` using the homomorphic encryption scheme `scheme`.
"""
function encrypt(msg::PlainText, scheme::BGVScheme)
    res = BGV(RLWE(scheme.oper.operQ.param.N, length(msg)), 0)
    encrypt_to!(res, msg, scheme)
    res
end

"""
    encrypt_to!(res::HE, msg::PlainText, scheme::HEScheme)

Encrypts a plaintext message `msg` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
function encrypt_to!(res::BGV, msg::PlainConst, scheme::BGVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert msg.isPQ[] "The input message should be in PQ."

    Qlen, auxQ = length(msg), msg.auxQ[]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = auxQ
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # Compute (tb + m, ta), which has the phase m + te mod Q.
    t = oper.ptxt_modulus.Q
    @inbounds for i = 1:Qlen
        _mul_to!(val.b.coeffs[i], t, val.b.coeffs[i], evalQ[i])
        _add_to!(val.b.coeffs[i], val.b.coeffs[i], msg.val[i], val.b.isntt[], evalQ[i])
        _mul_to!(val.a.coeffs[i], t, val.a.coeffs[i], evalQ[i])
    end

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i][1] == Qlen && Qatlevel[i][2] == auxQ
            res.level[] = i - 1
            return
        end
    end
    @error "There is no such level."
end

function encrypt_to!(res::BGV, msg::PlainPoly, scheme::BGVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert !msg.isPQ[] "The input message should not be in PQ."

    Qlen, auxQ = length(msg), msg.auxQ[]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = auxQ
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # NTT Buffer.
    buff_ntt = oper.tensor_buff[1].coeffs[1]

    # Compute (tb + m, ta), which has the phase m + te mod Q.
    t = oper.ptxt_modulus.Q
    @inbounds for i = 1:Qlen
        _mul_to!(res.val.b.coeffs[i], t, res.val.b.coeffs[i], evalQ[i])
        _mul_to!(res.val.a.coeffs[i], t, res.val.a.coeffs[i], evalQ[i])

        if !msg.val.isntt[]
            _ntt_to!(buff_ntt, msg.val.coeffs[i], evalQ[i])
            _add_to!(res.val.b.coeffs[i], res.val.b.coeffs[i], buff_ntt, evalQ[i])
        else
            _add_to!(res.val.b.coeffs[i], res.val.b.coeffs[i], msg.val.coeffs[i], evalQ[i])
        end
    end

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i][1] == Qlen && Qatlevel[i][2] == auxQ
            res.level[] = i - 1
            return
        end
    end
    @error "There is no such level."
end

"""
    decrypt(ct::HECiphertext, scheme::HEScheme)

Decrypts a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
function decrypt(ct::BGV, scheme::BGVScheme)
    res = PlainPoly(scheme.oper.operQ.param.N, length(ct.val))
    decrypt_to!(res, ct, scheme)
    res
end

"""
    decrypt_to!(res::PlainPoly, ct::HE, scheme::HEScheme)

Decrypts a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
function decrypt_to!(res::PlainPoly, ct::BGV, scheme::BGVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    resize!(res, length(ct.val))
    phase_to!(res, ct.val, entor)

    level = ct.level[]
    Qlen, auxQ = oper.Qatlevel[level+1]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)

    be = BasisExtender(evalQ.moduli, [oper.ptxt_modulus])
    @views basis_extend!(res.val.coeffs[1:1], res.val.coeffs, be)

    for i = Qlen:-1:1
        _Bred_to!(res.val.coeffs[i], res.val.coeffs[1], evalQ[i])
    end
end

"""
    error(ct::HECiphertext, scheme::HEScheme)

Computes the error of a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
function error(ct::BGV, scheme::BGVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    level = ct.level[]
    Qlen, auxQ = oper.Qatlevel[level+1]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)
    buff = oper.tensor_buff[2].coeffs[1:1]

    # Compute the phase.
    res = PlainPoly(scheme.oper.tensor_buff[1][1:Qlen])
    phase_to!(res, ct.val, entor)

    # Compute the error.
    be = BasisExtender(evalQ.moduli, [oper.ptxt_modulus])
    @views basis_extend!(buff, res.val.coeffs, be)
    buff = oper.tensor_buff[2].coeffs[1]
    for i = Qlen:-1:1
        tinvQi = invmod(oper.ptxt_modulus.Q, evalQ[i].Q)
        _sub_to!(res.val.coeffs[i], res.val.coeffs[i], buff, evalQ[i])
        _mul_to!(res.val.coeffs[i], tinvQi, res.val.coeffs[i], evalQ[i])
    end

    to_big(res.val, evalQ)
end

#================================================================================================================================#

"""
    BFVScheme(sketch::BFVParamSketch)
    BFVScheme(param::BFVParameters)
    BFVScheme(oper::BFVOperator)

BFVScheme is a struct for the BFV homomorphic encryption scheme. It contains the following fields: oper, entor, rlk, and atk.
- oper: BFVOperator
- entor: Union{Missing,Encryptor}
- rlk: Union{Missing,RLEV}
- atk: Dict{Int64,RLEV}
"""
mutable struct BFVScheme
    param::BFVParameters
    oper::BFVOperator
    entor::Union{Missing,Encryptor}
    rlk::Union{Missing,RLEV}
    atk::Dict{Int64,RLEV}

    BFVScheme(sketch::BFVParamSketch) = BFVScheme(BFVParameters(sketch))

    function BFVScheme(param::BFVParameters)
        oper = BFVOperator(param)

        new(param, oper, missing, missing, Dict{Int64,RLEV}())
    end
end

function encrypt(msg::PlainText, scheme::BFVScheme)
    res = BFV(RLWE(scheme.oper.operQ.param.N, length(msg)), 0)
    encrypt_to!(res, msg, scheme)
    res
end

function encrypt_to!(res::BFV, msg::PlainConst, scheme::BFVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert msg.isPQ[] "The input message should be in PQ."

    Qlen = length(msg)
    evalQ = geteval_at(Qlen, oper.operQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = 0
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # Compute (b + Δ⋅m, a) for Δ = ⌊Q/t⌉
    Δ = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q)
    @inbounds for i = 1:Qlen
        tmp = _mul(msg.val[i], (Δ % evalQ[i].Q) % UInt64, evalQ[i])
        _add_to!(val.b.coeffs[i], val.b.coeffs[i], tmp, val.b.isntt[], evalQ[i])
    end

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i] == Qlen
            res.level[] = i - 1
            return
        end
    end
    @error "There is no such level."
end

function encrypt_to!(res::BFV, msg::PlainPoly, scheme::BFVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert !msg.isPQ[] "The input message should not be in PQ."

    Qlen = length(msg)
    evalQ = geteval_at(Qlen, oper.operQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = 0
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # NTT Buffer.
    buff_ntt = oper.tensor_buff[1].coeffs[1]

    # Compute (b + Δ⋅m, a) for Δ = ⌊Q/t⌉.
    Δ = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q)
    @inbounds for i = 1:Qlen
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        if !msg.val.isntt[]
            _ntt_to!(buff_ntt, msg.val.coeffs[i], evalQ[i])
            _muladd_to!(res.val.b.coeffs[i], Δi, buff_ntt, evalQ[i])
        else
            _muladd_to!(res.val.b.coeffs[i], Δi, msg.val.coeffs[i], evalQ[i])
        end
    end

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i] == Qlen
            res.level[] = i - 1
            return
        end
    end
    @error "There is no such level."
end

function decrypt(ct::BFV, scheme::BFVScheme)
    res = PlainPoly(scheme.oper.operQ.param.N, length(ct.val))
    decrypt_to!(res, ct, scheme)
    res
end

function decrypt_to!(res::PlainPoly, ct::BFV, scheme::BFVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    resize!(res, length(ct.val))
    phase_to!(res, ct.val, entor)

    level = ct.level[]
    Qlen = oper.Qatlevel[level+1]
    evalQ = geteval_at(Qlen, oper.operQ)

    ss = SimpleScaler(evalQ.moduli, [oper.ptxt_modulus])
    @views simple_scale!(res.val.coeffs[1:1], res.val.coeffs, ss)

    for i = Qlen:-1:1
        _Bred_to!(res.val.coeffs[i], res.val.coeffs[1], evalQ[i])
    end
end

function error(ct::BFV, scheme::BFVScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    level = ct.level[]
    Qlen = oper.Qatlevel[level+1]
    evalQ = geteval_at(Qlen, oper.operQ)
    buff = oper.tensor_buff[2].coeffs[1:1]

    # Compute the phase.
    res = PlainPoly(scheme.oper.tensor_buff[1][1:Qlen])
    phase_to!(res, ct.val, entor)

    # Compute the error.
    ss = SimpleScaler(evalQ.moduli, [oper.ptxt_modulus])
    simple_scale!(buff, res.val.coeffs, ss)
    buff = oper.tensor_buff[2].coeffs[1]
    Δ = round(BigInt, prod(evalQ.moduli) // oper.ptxt_modulus.Q)
    for i = Qlen:-1:1
        Δi = (Δ % evalQ[i].Q.Q) % UInt64
        _mulsub_to!(res.val.coeffs[i], Δi, buff, evalQ[i])
    end

    to_big(res.val, evalQ)
end

#================================================================================================================================#

"""
    CKKSScheme(sketch::CKKSParamSketch)
    CKKSScheme(param::CKKSParameters)
    CKKSScheme(oper::CKKSOperator)

CKKSScheme is a struct for the CKKS homomorphic encryption scheme. It contains the following fields: oper, entor, rlk, and atk.
- oper: CKKSOperator
- entor: Union{Missing,Encryptor}
- rlk: Union{Missing,RLEV}
- atk: Dict{Int64,RLEV}
"""
mutable struct CKKSScheme
    oper::CKKSOperator
    entor::Union{Missing,Encryptor}
    rlk::Union{Missing,RLEV}
    atk::Dict{Int64,RLEV}

    CKKSScheme(sketch::CKKSParamSketch) = CKKSScheme(CKKSParameters(sketch))

    function CKKSScheme(param::CKKSParameters)
        oper = CKKSOperator(param)

        new(oper, missing, missing, Dict{Int64,RLEV}())
    end

    CKKSScheme(oper::CKKSOperator) = new(oper, missing, missing, Dict{Int64,RLEV}())
end

function encrypt(msg::PlainText, scheme::CKKSScheme)
    res = CKKS(RLWE(scheme.oper.operQ.param.N, length(msg)), 0)
    encrypt_to!(res, msg, scheme)
    res
end

function encrypt_to!(res::CKKS, msg::PlainConst, scheme::CKKSScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert msg.isPQ[] "The input message should be in PQ."

    Qlen, auxQ = length(msg), msg.auxQ[]
    evalQ = geteval_at(Qlen, oper.operQ, auxQ=auxQ)

    # Generate an RLWE sample, i.e., b + as = e mod Q.
    resize!(res.val, Qlen)
    res.val.auxQ[] = auxQ
    res.val.isPQ[] = false
    rlwe_sample_to!(res.val, entor)

    # Compute (b + m, a), which has the phase m + e mod Q.
    @inbounds for i = 1:Qlen
        _add_to!(val.b.coeffs[i], val.b.coeffs[i], msg.val[i], val.b.isntt[], evalQ[i])
    end

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i][1] == Qlen && Qatlevel[i][2] == auxQ
            res.level[] = i - 1
            return
        end
    end
    @error "There is no such level."
end

function encrypt_to!(res::CKKS, msg::PlainPoly, scheme::CKKSScheme)
    oper, entor = scheme.oper, scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert !msg.isPQ[] "The input message should not be in PQ."

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
            _ntt_to!(buff_ntt, msg.val.coeffs[i], evalQ[i])
            _add_to!(res.val.b.coeffs[i], res.val.b.coeffs[i], buff_ntt, evalQ[i])
        else
            _add_to!(res.val.b.coeffs[i], res.val.b.coeffs[i], msg.val.coeffs[i], evalQ[i])
        end
    end

    # Set level.
    Qatlevel = oper.Qatlevel
    for i = reverse(eachindex(Qatlevel))
        if Qatlevel[i][1] == Qlen && Qatlevel[i][2] == auxQ
            res.level[] = i - 1
            return
        end
    end
    @error "There is no such level."
end

function decrypt(ct::CKKS, scheme::CKKSScheme)
    res = PlainPoly(scheme.oper.operQ.param.N, length(ct.val))
    decrypt_to!(res, ct, scheme)
    res
end

function decrypt_to!(res::PlainPoly, ct::CKKS, scheme::CKKSScheme)
    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    resize!(res, length(ct.val))
    phase_to!(res, ct.val, entor)
end

#================================================================================================================================#

const IntScheme = Union{BFVScheme,BGVScheme}
const HEScheme = Union{BFVScheme,BGVScheme,CKKSScheme}
const HEOperator = Union{BFVOperator,BGVOperator,CKKSOperator}
const HECiphertext = Union{BFV,BGV,CKKS}

export BGVScheme, BFVScheme, CKKSScheme
export encrypt, encrypt_to!, decrypt, decrypt_to!, error
export IntScheme, HEScheme, HEOperator, HECiphertext

#================================================================================#

"""
    set_encryptor!(key::RLWEkey, scheme::HEScheme; σ::Real=3.2)
    set_encryptor!(pk::RLWE, scheme::HEScheme; σ::Real=3.2, τ::Real=3.2)

Sets an encryptor for the homomorphic encryption scheme `scheme` using the secret key `key` or the public key `pk`.
For SKEncryptor, σ denotes the standard deviation of the discrete Gaussian distribution of the noise.
For PKEncryptor, σ and τ denote the standard deviation for the discrete Gaussian distribution multiplied to the public key, 
and the standard deviatino for the rounded Gaussian distribution for the additional noise, respectively.
"""
set_encryptor!(key::RLWEkey, scheme::HEScheme; σ::Real=3.2) = begin
    operQ = scheme.oper.operQ
    entor = SKEncryptor(key, σ, operQ)
    scheme.entor = entor
end
set_encryptor!(pk::RLWE, scheme::HEScheme; σ::Real=3.2, τ::Real=3.2) = begin
    operQ = scheme.oper.operQ
    entor = PKEncryptor(pk, σ, τ, operQ)
    scheme.entor = entor
end

"""
    public_keygen(scheme::HEScheme)

Generates a public key for the homomorphic encryption scheme `scheme`.
The scheme should have a secret key encryptor generated.
"""
public_keygen(scheme::HEScheme) = begin
    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    public_keygen(entor)
end

"""
    relin_keygen(scheme::HEScheme)

Generates a relinearization key for the homomorphic encryption scheme `scheme`.
The scheme should have a secret key encryptor generated.
"""
relin_keygen(scheme::HEScheme) = begin
    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    relin_keygen(entor)
end
"""
    set_relinkey!(rlk::RLEV, scheme::HEScheme)

Sets a relinearization key `rlk` for the homomorphic encryption scheme `scheme`.
"""
set_relinkey!(rlk::RLEV, scheme::HEScheme) = (scheme.rlk = rlk)

"""
    rotate_keygen(idx::NTuple{N,Int64}, scheme::HEScheme)
    rotate_keygen(idxset::Vector{NTuple{N,Int64}}, scheme::HEScheme)

Generates a rotation key for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
The scheme should have a secret key encryptor generated.
"""
rotate_keygen(idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_keygen(autidx, entor)
end
rotate_keygen(idxset::Vector{NTuple{N,Int64}}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    rtks = Vector{RLEV}(undef, length(idxset))
    for i = eachindex(idxset)
        idx = idxset[i]
        autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
        rtks[i] = automorphism_keygen(autidx, entor)
    end

    rtks
end

"""
    set_rotate_key!(idx::NTuple{N,Int64}, rtk::RLEV, scheme::HEScheme)
    set_rotate_key!(idxset::Vector{NTuple{N,Int64}}, rtkset::Vector{RLEV}, scheme::HEScheme)

Sets a rotation key `rtk` for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
"""
set_rotate_key!(idx::NTuple{N,Int64}, rtk::RLEV, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)

    scheme.atk[autidx] = rtk
end
set_rotate_key!(idxset::Vector{NTuple{N,Int64}}, rtkset::Vector{RLEV}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    for (idx, rtk) ∈ zip(idxset, rtkset)
        autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
        scheme.atk[autidx] = rtk
    end
end

"""
    automorphism_keygen(idx::Int64, scheme::HEScheme)
    automorphism_keygen(idxset::Vector{Int64}, scheme::HEScheme)

Generates an automorphism key for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
The scheme should have a secret key encryptor generated.
"""
automorphism_keygen(idx::Int64, scheme::HEScheme) = begin
    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    automorphism_keygen(idx, entor)
end
automorphism_keygen(idxset::Vector{Int64}, scheme::HEScheme) = begin
    entor = scheme.entor
    @assert !ismissing(entor) "The encryptor is not generated."
    @assert typeof(entor) == SKEncryptor "The Encryptor should be a secret key encryptor."

    [automorphism_keygen(idx, entor) for idx ∈ idxset]
end

"""
    set_automorphism_key!(idx::Int64, atk::RLEV, scheme::HEScheme)
    set_automorphism_key!(idxset::Vector{Int64}, atkset::Vector{RLEV}, scheme::HEScheme)

Sets an automorphism key `atk` for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
"""
set_automorphism_key!(idx::Int64, atk::RLEV, scheme::HEScheme) = (scheme.atk[idx] = atk)
set_automorphism_key!(idxset::Vector{Int64}, atkset::Vector{RLEV}, scheme::HEScheme) = begin
    for (idx, atk) ∈ zip(idxset, atkset)
        scheme.atk[idx] = atk
    end
end

"""
    encode(msg::Union{Int64,UInt64}, scheme::IntScheme; level::Integer=typemax(Int64), isPQ::Bool=false)
    encode(msg::Vector{UInt64}, scheme::IntScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true)
    encode(msg::AbstractFloat, scheme::CKKSScheme; level::Integer=typemax(Int64), isPQ::Bool=false)
    encode(msg::Vector{<:Union{AbstractFloat,Complex{<:AbstractFloat}}}, scheme::CKKSScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true)

Encodes a message `msg` using the homomorphic encryption scheme `scheme`.
"""
encode(msg::UInt64, scheme::IntScheme; level::Integer=typemax(Int64), isPQ::Bool=false) = encode(msg, scheme.oper, level=level, isPQ=isPQ)
encode(msg::Vector{UInt64}, scheme::IntScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true) = encode(msg, scheme.oper, level=level, isPQ=isPQ, ispacking=ispacking)
encode(msg::AbstractFloat, scheme::CKKSScheme; level::Integer=typemax(Int64), isPQ::Bool=false) = encode(msg, scheme.oper, level=level, isPQ=isPQ)
encode(msg::Vector{<:Union{AbstractFloat,Complex{<:AbstractFloat}}}, scheme::CKKSScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true) = encode(msg, scheme.oper, level=level, isPQ=isPQ, ispacking=ispacking)

"""
    decode(pt::PlainPoly, scheme::HEScheme; ispacking::Bool=true)
    decode(pt::PlainConst, scheme::HEScheme)

Decodes a plaintext `pt` using the homomorphic encryption scheme `scheme`.
"""
decode(pt::PlainConst, scheme::HEScheme) = decode(pt, scheme.oper)
decode(pt::PlainPoly, scheme::HEScheme; ispacking::Bool=true) = decode(pt, scheme.oper, ispacking=ispacking)

"""
    change_level(pt::PlainText, targetlvl::Integer, scheme::HEScheme)
    
Changes the level of a plaintext `pt` to the target level `targetlvl` using the homomorphic encryption scheme `scheme`.
"""
change_level(pt::PlainText, targetlvl::Integer, scheme::HEScheme) = change_level(pt, targetlvl, scheme.oper)
"""
    change_level_to!(pt::PlainText, targetlvl::Integer, scheme::HEScheme)

Changes the level of a plaintext `pt` to the target level `targetlvl` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
change_level_to!(res::PlainText, pt::PlainText, targetlvl::Integer, scheme::HEScheme) = change_level_to!(res, pt, targetlvl, scheme.oper)

"""
    drop_level_to!(res::HECiphertext, ct::HECiphertext, targetlvl::Integer, scheme::HEScheme)

Drops the level of a ciphertext `ct` to the target level `targetlvl` and stores the result in `res`.
"""
drop_level_to!(res::HECiphertext, ct::HECiphertext, targetlvl::Integer, scheme::HEScheme) = drop_level_to!(res, ct, targetlvl, scheme.oper)

"""
    rescale(ct::HECiphertext, scheme::HEScheme)

Rescales a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
rescale(ct::HECiphertext, scheme::HEScheme) = begin
    res = similar(ct)
    rescale_to!(res, ct, scheme)
    res
end

"""
    rescale_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme)

Rescales a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
rescale_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme) = rescale_to!(res, ct, scheme.oper)

"""
    neg(ct::HECiphertext, scheme::HEScheme)

Negates a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
neg(ct::HECiphertext, scheme::HEScheme) = begin
    res = similar(ct)
    neg_to!(res, ct, scheme)
    res
end
"""
    neg_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme)

Negates a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
neg_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme) = neg_to!(res, ct, scheme.oper)

"""
    add(x::HECiphertext, y::PlainPoly, scheme::HEScheme)
    add(x::PlainPoly, y::HECiphertext, scheme::HEScheme)
    add(x::HECiphertext, y::PlainConst, scheme::HEScheme)
    add(x::PlainConst, y::HECiphertext, scheme::HEScheme)
    add(x::HECiphertext, y::HECiphertext, scheme::HEScheme)

Adds two operands `x` and `y` using the homomorphic encryption scheme `scheme`.
"""
add(x::HECiphertext, y::PlainPoly, scheme::HEScheme) = begin
    res = similar(x)
    add_to!(res, x, y, scheme)
    res
end
add_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme) = add_to!(res, x, y, scheme.oper)
add(x::PlainPoly, y::HECiphertext, scheme::HEScheme) = add(y, x, scheme)
add_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme) = add_to!(res, y, x, scheme)

add(x::HECiphertext, y::PlainConst, scheme::HEScheme) = begin
    res = similar(x)
    add_to!(res, x, y, scheme)
    res
end
add_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme) = add_to!(res, x, y, scheme.oper)
add(x::PlainConst, y::HECiphertext, scheme::HEScheme) = add(y, x, scheme)
add_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme) = add_to!(res, y, x, scheme)

add(x::HECiphertext, y::HECiphertext, scheme::HEScheme) = begin
    res = similar(x)
    add_to!(res, x, y, scheme)
    res
end
add_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme) = add_to!(res, x, y, scheme.oper)

"""
    sub(x::HECiphertext, y::PlainPoly, scheme::HEScheme)
    sub(x::PlainPoly, y::HECiphertext, scheme::HEScheme)
    sub(x::HECiphertext, y::PlainConst, scheme::HEScheme)
    sub(x::PlainConst, y::HECiphertext, scheme::HEScheme)
    sub(x::HECiphertext, y::HECiphertext, scheme::HEScheme)

Subtracts two operands `x` and `y` using the homomorphic encryption scheme `scheme`.
"""
sub(x::HECiphertext, y::PlainPoly, scheme::HEScheme) = begin
    res = similar(x)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme) = sub_to!(res, x, y, scheme.oper)
sub(x::PlainPoly, y::HECiphertext, scheme::HEScheme) = begin
    res = similar(y)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme) = sub_to!(res, x, y, scheme.oper)

sub(x::HECiphertext, y::PlainConst, scheme::HEScheme) = begin
    res = similar(x)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme) = sub_to!(res, x, y, scheme.oper)
sub(x::PlainConst, y::HECiphertext, scheme::HEScheme) = begin
    res = similar(y)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme) = sub_to!(res, x, y, scheme.oper)

sub(x::HECiphertext, y::HECiphertext, scheme::HEScheme) = begin
    res = similar(x)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme) = sub_to!(res, x, y, scheme.oper)

# MUL
"""
    mul(x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false)
    mul(x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)
    mul(x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false)
    mul(x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)
    mul(x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)

Multiplies two operands `x` and `y` using the homomorphic encryption scheme `scheme`.
"""
mul(x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, scheme, islazy=islazy)
    res
end
"""
    mul_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false)
    mul_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)
    mul_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false)
    mul_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)
    mul_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)

Multiplies two operands `x` and `y` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
mul_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false) = mul_to!(res, x, y, scheme.oper, islazy=islazy)
mul(x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false) = mul(y, x, scheme, islazy=islazy)
mul_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false) = mul_to!(res, y, x, scheme, islazy=islazy)

mul(x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, scheme, islazy=islazy)
    res
end
mul_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false) = mul_to!(res, x, y, scheme.oper, islazy=islazy)
mul(x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false) = mul(y, x, scheme, islazy=islazy)
mul_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false) = mul_to!(res, y, x, scheme, islazy=islazy)

mul(x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false) = begin
    res = similar(x)
    mul_to!(res, x, y, scheme, islazy=islazy)
    res
end
mul_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false) = begin
    @assert !ismissing(scheme.rlk) "The relinearization key is not generated."
    mul_to!(res, x, y, scheme.rlk, scheme.oper, islazy=islazy)
end

"""
    keyswitch(ct::HECiphertext, ksk::RLEV, scheme::HEScheme)

Performs a keyswitch operation on a ciphertext `ct` using the key-switching key `ksk`.
"""
keyswitch(ct::HECiphertext, ksk::RLEV, scheme::HEScheme) = begin
    res = similar(ct)
    keyswitch_to!(res, ct, ksk, scheme)
    res
end
"""
    keyswitch_to!(res::HECiphertext, ct::HECiphertext, ksk::RLEV, scheme::HEScheme)

Performs a keyswitch operation on a ciphertext `ct` using the key-switching key `ksk` and stores the result in `res`.
"""
keyswitch_to!(res::HECiphertext, ct::HECiphertext, ksk::RLEV, scheme::HEScheme) = keyswitch_to!(res, ct, ksk, scheme.oper)
"""
    rotate(ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
rotate(ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    res = similar(ct)
    rotate_to!(res, ct, idx, scheme)
    res
end
"""
    rotate_to!(res::HECiphertext, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
rotate_to!(res::HECiphertext, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, ct, autidx, scheme)
end
"""
    hoisted_rotate(adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a hoisted rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
hoisted_rotate(adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    res = similar(ct)
    hoisted_rotate_to!(res, adec, ct, idx, scheme)
    res
end
"""
    hoisted_rotate_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a hoisted rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
hoisted_rotate_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    @assert !ismissing(packer) "Rotation operation cannot be defined without SIMD packing."

    cube, cubegen = packer.cube, packer.cubegen
    @assert length(cube) == N "The number of indices should match the number of dimensions."

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, ct, autidx, scheme)
end
"""
    automorphism(ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies an automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
automorphism(ct::HECiphertext, idx::Int64, scheme::HEScheme) = begin
    res = similar(ct)
    automorphism_to!(res, ct, idx, scheme)
    res
end
"""
    automorphism_to!(res::HECiphertext, ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies an automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
automorphism_to!(res::HECiphertext, ct::HECiphertext, idx::Int64, scheme::HEScheme) = begin
    @assert haskey(scheme.atk, idx) "The rotation key is not generated."
    automorphism_to!(res, ct, idx, scheme.atk[idx], scheme.oper)
end
"""
    hoisted_automorphism(adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies a hoisted automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
hoisted_automorphism(adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme) = begin
    res = similar(ct)
    hoisted_automorphism_to!(res, adec, ct, idx, scheme)
    res
end
"""
    hoisted_automorphism_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies a hoisted automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
hoisted_automorphism_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme) = begin
    @assert haskey(scheme.atk, idx) "The rotation key is not generated."
    hoisted_automorphism_to!(res, adec, ct, idx, scheme.atk[idx], scheme.oper)
end

export set_encryptor!, public_keygen, relin_keygen, set_relinkey!, rotate_keygen, set_rotate_key!, rotate_keygen, set_rotate_key!, automorphism_keygen, set_automorphism_key!, automorphism_keygen, set_automorphism_key!