"""
    HEScheme

Abstract type for FHE schemes.
"""
abstract type HEScheme end

"""
    HEParamSketch

Abstract type for FHE parameter sketches.
"""
abstract type HEParamSketch end

"""
    HEParameters

Abstract type for FHE parameters.
"""
abstract type HEParameters end

"""
    HEOperator

Abstract type for FHE operators.
"""
abstract type HEOperator end

"""
    HECiphertext

Abstract type for FHE ciphertexts.
"""
abstract type HECiphertext end

#================================================================================#

"""
    set_encryptor!(key::RLWEkey, scheme::HEScheme; σ::Real=3.2)
    set_encryptor!(pk::RLWE, scheme::HEScheme; σ::Real=3.2, τ::Real=3.2)

Sets an encryptor for the homomorphic encryption scheme `scheme` using the secret key `key` or the public key `pk`.
For SKEncryptor, σ denotes the standard deviation of the discrete Gaussian distribution of the noise.
For PKEncryptor, σ and τ denote the standard deviation for the discrete Gaussian distribution multiplied to the public key, 
and the standard deviatino for the rounded Gaussian distribution for the additional noise, respectively.
"""
set_encryptor!(key::RLWEkey, scheme::HEScheme; σ::Real=3.2)::Nothing = begin
    operQ = scheme.oper.operQ
    entor = SKEncryptor(key, σ, operQ)
    scheme.entor = entor
    return nothing
end
set_encryptor!(pk::RLWE, scheme::HEScheme; σ::Real=3.2, τ::Real=3.2)::Nothing = begin
    operQ = scheme.oper.operQ
    entor = PKEncryptor(pk, σ, τ, operQ)
    scheme.entor = entor
    return nothing
end

"""
    encrypt(msg::PlainText, scheme::HEScheme)

Encrypts a plaintext message `msg` using the homomorphic encryption scheme `scheme`.
"""
encrypt(msg::PlainText, scheme::HEScheme)::HECiphertext = begin
    entor, oper = scheme.entor, scheme.oper
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if msg.isPQ[]
        throw(DomainError("The input message should not be in PQ."))
    end

    encrypt(msg, entor, oper)
end

"""
    encrypt_to!(res::HECiphertext, msg::PlainText, scheme::HEScheme)

Encrypts a plaintext message `msg` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
encrypt_to!(res::HECiphertext, msg::PlainText, scheme::HEScheme)::Nothing = begin
    entor, oper = scheme.entor, scheme.oper
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if msg.isPQ[]
        throw(DomainError("The input message should not be in PQ."))
    end

    encrypt_to!(res, msg, entor, oper)

    return nothing
end

"""
    decrypt(ct::HECiphertext, scheme::HEScheme)

Decrypts a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
decrypt(ct::HECiphertext, scheme::HEScheme)::PlainText = begin
    entor, oper = scheme.entor, scheme.oper
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end
    decrypt(ct, entor, oper)
end

"""
    decrypt_to!(res::PlainText, ct::HECiphertext, scheme::HEScheme)

Decrypts a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
decrypt_to!(res::PlainText, ct::HECiphertext, scheme::HEScheme)::Nothing = begin
    entor, oper = scheme.entor, scheme.oper
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end
    decrypt_to!(res, ct, entor, oper)

    return nothing
end

"""
    get_error(ct::HECiphertext, scheme::HEScheme)

Computes the error of a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
get_error(ct::HECiphertext, scheme::HEScheme)::Vector{BigInt} = begin
    entor, oper = scheme.entor, scheme.oper
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end
    get_error(ct, entor, oper)
end

"""
    public_keygen(scheme::HEScheme)

Generates a public key for the homomorphic encryption scheme `scheme`.
The scheme should have a secret key encryptor generated.
"""
public_keygen(scheme::HEScheme)::RLWE = begin
    entor = scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    public_keygen(entor)
end

"""
    relin_keygen(scheme::HEScheme)

Generates a relinearization key for the homomorphic encryption scheme `scheme`.
The scheme should have a secret key encryptor generated.
"""
relin_keygen(scheme::HEScheme)::RLEV = begin
    entor = scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    relin_keygen(entor)
end
"""
    set_relinkey!(rlk::RLEV, scheme::HEScheme)

Sets a relinearization key `rlk` for the homomorphic encryption scheme `scheme`.
"""
set_relinkey!(rlk::RLEV, scheme::HEScheme)::Nothing = begin
    (scheme.rlk = rlk)
    return nothing
end

"""
    rotate_keygen(idx::NTuple{N,Int64}, scheme::HEScheme)
    rotate_keygen(idxset::Vector{NTuple{N,Int64}}, scheme::HEScheme)

Generates a rotation key for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
The scheme should have a secret key encryptor generated.
"""
rotate_keygen(idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    entor = scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_keygen(autidx, entor)
end::RLEV
rotate_keygen(idxset::Vector{NTuple{N,Int64}}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    entor = scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    rtks = Vector{RLEV}(undef, length(idxset))
    for i = eachindex(idxset)
        idx = idxset[i]
        autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
        rtks[i] = automorphism_keygen(autidx, entor)
    end

    rtks
end::Vector{RLEV}

"""
    set_rotate_key!(idx::NTuple{N,Int64}, rtk::RLEV, scheme::HEScheme)
    set_rotate_key!(idxset::Vector{NTuple{N,Int64}}, rtkset::Vector{RLEV}, scheme::HEScheme)

Sets a rotation key `rtk` for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
"""
set_rotate_key!(idx::NTuple{N,Int64}, rtk::RLEV, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)

    scheme.atk[autidx] = rtk

    return nothing
end::Nothing
set_rotate_key!(idxset::Vector{NTuple{N,Int64}}, rtkset::Vector{RLEV}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    for (idx, rtk) ∈ zip(idxset, rtkset)
        autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
        scheme.atk[autidx] = rtk
    end

    return nothing
end::Nothing

"""
    automorphism_keygen(idx::Int64, scheme::HEScheme)
    automorphism_keygen(idxset::Vector{Int64}, scheme::HEScheme)

Generates an automorphism key for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
The scheme should have a secret key encryptor generated.
"""
automorphism_keygen(idx::Int64, scheme::HEScheme)::RLEV = begin
    entor = scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    automorphism_keygen(idx, entor)
end
automorphism_keygen(idxset::Vector{Int64}, scheme::HEScheme)::Vector{RLEV} = begin
    entor = scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    [automorphism_keygen(idx, entor) for idx ∈ idxset]
end

"""
    set_automorphism_key!(idx::Int64, atk::RLEV, scheme::HEScheme)
    set_automorphism_key!(idxset::Vector{Int64}, atkset::Vector{RLEV}, scheme::HEScheme)

Sets an automorphism key `atk` for the homomorphic encryption scheme `scheme` with the index `idx`, or index set `idxset`.
"""
set_automorphism_key!(idx::Int64, atk::RLEV, scheme::HEScheme)::Nothing = begin
    (scheme.atk[idx] = atk)
    return nothing
end
set_automorphism_key!(idxset::Vector{Int64}, atkset::Vector{RLEV}, scheme::HEScheme)::Nothing = begin
    for (idx, atk) ∈ zip(idxset, atkset)
        scheme.atk[idx] = atk
    end

    return nothing
end

"""
    encode(msg::Number, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, scaling_factor::Real=0.0)
    encode(msg::Vector{<:Number}, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true, scaling_factor::Real=0.0)

Encodes a message `msg` using the homomorphic encryption scheme `scheme`.
The scaling factor is only used for approximate homomorphic schemes.
"""
encode(msg::Number, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, scaling_factor::Real=0.0)::PlainConst =
    encode(msg, scheme.oper, level=level, isPQ=isPQ, scaling_factor=scaling_factor)
encode(msg::Vector{<:Number}, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true, scaling_factor::Real=0.0)::PlainPoly =
    encode(msg, scheme.oper, level=level, isPQ=isPQ, ispacking=ispacking, scaling_factor=scaling_factor)

"""
    encode_to!(res::PlainConst, msg::Number, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, scaling_factor::Real=0.0)
    encode_to!(res::PlainPoly, msg::Vector{<:Number}, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true, scaling_factor::Real=0.0)

Encodes a message `msg` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
The scaling factor is only used for approximate homomorphic schemes.
"""
encode_to!(res::PlainConst, msg::Number, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, scaling_factor::Real=0.0)::Nothing =
    encode_to!(res, msg, scheme.oper, level=level, isPQ=isPQ, scaling_factor=scaling_factor)
encode_to!(res::PlainPoly, msg::Vector{<:Number}, scheme::HEScheme; level::Integer=typemax(Int64), isPQ::Bool=false, ispacking::Bool=true, scaling_factor::Real=0.0)::Nothing =
    encode_to!(res, msg, scheme.oper, level=level, isPQ=isPQ, ispacking=ispacking, scaling_factor=scaling_factor)

"""
    decode(pt::PlainPoly, scheme::HEScheme; ispacking::Bool=true)
    decode(pt::PlainConst, scheme::HEScheme)

Decodes a plaintext `pt` using the homomorphic encryption scheme `scheme`.
"""
decode(pt::PlainConst, scheme::HEScheme)::Number = decode(pt, scheme.oper)
decode(pt::PlainPoly, scheme::HEScheme; ispacking::Bool=true)::Vector{<:Number} = decode(pt, scheme.oper, ispacking=ispacking)

"""
    decode_to!(res::Vector{<:Number}, pt::PlainPoly, scheme::HEScheme; ispacking::Bool=true)

Decodes a plaintext `pt` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
decode_to!(res::Vector{<:Number}, pt::PlainPoly, scheme::HEScheme; ispacking::Bool=true)::Nothing = decode_to!(res, pt, scheme.oper, ispacking=ispacking)

"""
    change_level(pt::PlainText, targetlvl::Integer, scheme::HEScheme)
    
Changes the level of a plaintext `pt` to the target level `targetlvl` using the homomorphic encryption scheme `scheme`.
"""
change_level(pt::PlainText, targetlvl::Integer, scheme::HEScheme)::PlainText = change_level(pt, targetlvl, scheme.oper)
"""
    change_level_to!(pt::PlainText, targetlvl::Integer, scheme::HEScheme)

Changes the level of a plaintext `pt` to the target level `targetlvl` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
change_level_to!(res::PlainText, pt::PlainText, targetlvl::Integer, scheme::HEScheme)::Nothing = begin
    change_level_to!(res, pt, targetlvl, scheme.oper)
    return nothing
end
"""
    drop_level_to!(res::HECiphertext, ct::HECiphertext, targetlvl::Integer, scheme::HEScheme)

Drops the level of a ciphertext `ct` to the target level `targetlvl` and stores the result in `res`.
"""
drop_level_to!(res::HECiphertext, ct::HECiphertext, targetlvl::Integer, scheme::HEScheme)::Nothing = begin
    drop_level_to!(res, ct, targetlvl, scheme.oper)
    return nothing
end

"""
    rescale(ct::HECiphertext, scheme::HEScheme)

Rescales a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
rescale(ct::HECiphertext, scheme::HEScheme)::HECiphertext = begin
    res = similar(ct)
    rescale_to!(res, ct, scheme)
    res
end

"""
    rescale_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme)

Rescales a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
rescale_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme)::Nothing = begin
    rescale_to!(res, ct, scheme.oper)
    return nothing
end

"""
    neg(ct::HECiphertext, scheme::HEScheme)

Negates a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
neg(ct::HECiphertext, scheme::HEScheme)::HECiphertext = begin
    res = similar(ct)
    neg_to!(res, ct, scheme)
    res
end
"""
    neg_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme)

Negates a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
neg_to!(res::HECiphertext, ct::HECiphertext, scheme::HEScheme)::Nothing = begin
    neg_to!(res, ct, scheme.oper)
    return nothing
end

"""
    add(x::HECiphertext, y::PlainPoly, scheme::HEScheme)
    add(x::PlainPoly, y::HECiphertext, scheme::HEScheme)
    add(x::HECiphertext, y::PlainConst, scheme::HEScheme)
    add(x::PlainConst, y::HECiphertext, scheme::HEScheme)
    add(x::HECiphertext, y::HECiphertext, scheme::HEScheme)

Adds two operands `x` and `y` using the homomorphic encryption scheme `scheme`.
"""
add(x::HECiphertext, y::PlainPoly, scheme::HEScheme)::HECiphertext = begin
    res = similar(x)
    add_to!(res, x, y, scheme)
    res
end
add_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme)::Nothing = begin
    add_to!(res, x, y, scheme.oper)
    return nothing
end
add(x::PlainPoly, y::HECiphertext, scheme::HEScheme)::HECiphertext = add(y, x, scheme)
add_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme)::Nothing = begin
    add_to!(res, y, x, scheme)
    return nothing
end

add(x::HECiphertext, y::PlainConst, scheme::HEScheme)::HECiphertext = begin
    res = similar(x)
    add_to!(res, x, y, scheme)
    res
end
add_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme)::Nothing = begin
    add_to!(res, x, y, scheme.oper)
    return nothing
end
add(x::PlainConst, y::HECiphertext, scheme::HEScheme)::HECiphertext = add(y, x, scheme)
add_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme)::Nothing = begin
    add_to!(res, y, x, scheme)
    return nothing
end

add(x::HECiphertext, y::HECiphertext, scheme::HEScheme)::HECiphertext = begin
    res = similar(x)
    add_to!(res, x, y, scheme)
    res
end
add_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme)::Nothing = begin
    add_to!(res, x, y, scheme.oper)
    return nothing
end

"""
    sub(x::HECiphertext, y::PlainPoly, scheme::HEScheme)
    sub(x::PlainPoly, y::HECiphertext, scheme::HEScheme)
    sub(x::HECiphertext, y::PlainConst, scheme::HEScheme)
    sub(x::PlainConst, y::HECiphertext, scheme::HEScheme)
    sub(x::HECiphertext, y::HECiphertext, scheme::HEScheme)

Subtracts two operands `x` and `y` using the homomorphic encryption scheme `scheme`.
"""
sub(x::HECiphertext, y::PlainPoly, scheme::HEScheme)::HECiphertext = begin
    res = similar(x)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme)::Nothing = begin
    sub_to!(res, x, y, scheme.oper)
    return nothing
end
sub(x::PlainPoly, y::HECiphertext, scheme::HEScheme)::HECiphertext = begin
    res = similar(y)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme)::Nothing = begin
    sub_to!(res, x, y, scheme.oper)
    return nothing
end

sub(x::HECiphertext, y::PlainConst, scheme::HEScheme)::HECiphertext = begin
    res = similar(x)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme)::Nothing = begin
    sub_to!(res, x, y, scheme.oper)
    return nothing
end
sub(x::PlainConst, y::HECiphertext, scheme::HEScheme)::HECiphertext = begin
    res = similar(y)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme)::Nothing = begin
    sub_to!(res, x, y, scheme.oper)
    return nothing
end

sub(x::HECiphertext, y::HECiphertext, scheme::HEScheme)::HECiphertext = begin
    res = similar(x)
    sub_to!(res, x, y, scheme)
    res
end
sub_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme)::Nothing = begin
    sub_to!(res, x, y, scheme.oper)
    return nothing
end

# MUL
"""
    mul(x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false)
    mul(x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)
    mul(x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false)
    mul(x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)
    mul(x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)

Multiplies two operands `x` and `y` using the homomorphic encryption scheme `scheme`.
"""
mul(x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false)::HECiphertext = begin
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
mul_to!(res::HECiphertext, x::HECiphertext, y::PlainPoly, scheme::HEScheme; islazy::Bool=false)::Nothing = begin
    mul_to!(res, x, y, scheme.oper, islazy=islazy)
    return nothing
end
mul(x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)::HECiphertext = begin
    res = similar(y)
    mul_to!(res, y, x, scheme, islazy=islazy)
    res
end
mul_to!(res::HECiphertext, x::PlainPoly, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)::Nothing = begin
    mul_to!(res, y, x, scheme.oper, islazy=islazy)
    return nothing
end

mul(x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false)::HECiphertext = begin
    res = similar(x)
    mul_to!(res, x, y, scheme, islazy=islazy)
    res
end
mul_to!(res::HECiphertext, x::HECiphertext, y::PlainConst, scheme::HEScheme; islazy::Bool=false)::Nothing = begin
    mul_to!(res, x, y, scheme.oper, islazy=islazy)
    return nothing
end
mul(x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)::HECiphertext = begin
    res = similar(y)
    mul_to!(res, y, x, scheme, islazy=islazy)
    res
end
mul_to!(res::HECiphertext, x::PlainConst, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)::Nothing = begin
    mul_to!(res, y, x, scheme.oper, islazy=islazy)
    return nothing
end

mul(x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)::HECiphertext = begin
    res = similar(x)
    mul_to!(res, x, y, scheme, islazy=islazy)
    res
end
mul_to!(res::HECiphertext, x::HECiphertext, y::HECiphertext, scheme::HEScheme; islazy::Bool=false)::Nothing = begin
    if ismissing(scheme.rlk)
        throw(ErrorException("The relinearization key is not generated."))
    end
    mul_to!(res, x, y, scheme.rlk, scheme.oper, islazy=islazy)
    return nothing
end

"""
    keyswitch(ct::HECiphertext, ksk::RLEV, scheme::HEScheme)

Performs a keyswitch operation on a ciphertext `ct` using the key-switching key `ksk`.
"""
keyswitch(ct::HECiphertext, ksk::RLEV, scheme::HEScheme)::HECiphertext = begin
    res = similar(ct)
    keyswitch_to!(res, ct, ksk, scheme)
    res
end
"""
    keyswitch_to!(res::HECiphertext, ct::HECiphertext, ksk::RLEV, scheme::HEScheme)

Performs a keyswitch operation on a ciphertext `ct` using the key-switching key `ksk` and stores the result in `res`.
"""
keyswitch_to!(res::HECiphertext, ct::HECiphertext, ksk::RLEV, scheme::HEScheme)::Nothing = begin
    keyswitch_to!(res, ct, ksk, scheme.oper)
    return nothing
end
"""
    rotate(ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
rotate(ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    res = similar(ct)
    rotate_to!(res, ct, idx, scheme)
    res
end::HECiphertext
"""
    rotate_to!(res::HECiphertext, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
rotate_to!(res::HECiphertext, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    automorphism_to!(res, ct, autidx, scheme)

    return nothing
end::Nothing
"""
    hoisted_rotate(adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a hoisted rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
hoisted_rotate(adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    res = similar(ct)
    hoisted_rotate_to!(res, adec, ct, idx, scheme)
    res
end::HECiphertext
"""
    hoisted_rotate_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N}

Applies a hoisted rotation to a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
hoisted_rotate_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::NTuple{N,Int64}, scheme::HEScheme) where {N} = begin
    packer = scheme.oper.packer
    if ismissing(packer)
        throw(ErrorException("Rotation operation cannot be defined without SIMD packing."))
    end

    cube, cubegen = packer.cube, packer.cubegen
    if length(cube) ≠ N
        throw(DimensionMismatch("The number of indices should match the number of dimensions."))
    end

    autidx = gen_power_modmul(cubegen, cube, idx, packer.m)
    hoisted_automorphism_to!(res, adec, ct, autidx, scheme)

    return nothing
end::Nothing
"""
    automorphism(ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies an automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
automorphism(ct::HECiphertext, idx::Int64, scheme::HEScheme)::HECiphertext = begin
    res = similar(ct)
    automorphism_to!(res, ct, idx, scheme)
    res
end
"""
    automorphism_to!(res::HECiphertext, ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies an automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
automorphism_to!(res::HECiphertext, ct::HECiphertext, idx::Int64, scheme::HEScheme)::Nothing = begin
    if !haskey(scheme.atk, idx)
        throw(ErrorException("The rotation key is not generated."))
    end
    automorphism_to!(res, ct, idx, scheme.atk[idx], scheme.oper)

    return nothing
end
"""
    hoisted_automorphism(adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies a hoisted automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme`.
"""
hoisted_automorphism(adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme)::HECiphertext = begin
    res = similar(ct)
    hoisted_automorphism_to!(res, adec, ct, idx, scheme)
    res
end
"""
    hoisted_automorphism_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme)

Applies a hoisted automorphism to a ciphertext `ct` using the homomorphic encryption scheme `scheme` and stores the result in `res`.
"""
hoisted_automorphism_to!(res::HECiphertext, adec::Tensor, ct::HECiphertext, idx::Int64, scheme::HEScheme)::Nothing = begin
    if !haskey(scheme.atk, idx)
        throw(ErrorException("The rotation key is not generated."))
    end
    hoisted_automorphism_to!(res, adec, ct, idx, scheme.atk[idx], scheme.oper)

    return nothing
end
