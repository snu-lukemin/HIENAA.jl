"""
    LeveledScheme

Internal interface for the leveled HE schemes.
"""
module LeveledScheme

import ..Math: add, add_to!, sub, sub_to!, neg, neg_to!
import ..Math: gen_power_modmul

import ..Ring: mul, mul_to!, automorphism, automorphism_to!, automorphism!, initialise!

import ..Rlwe: PlainConst, PlainPoly, PlainText, RLWE, Tensor, RLEV, decompose
import ..Rlwe: scale_to!, keyswitch, keyswitch_to!, hoisted_keyswitch, hoisted_keyswitch_to!, hoisted_automorphism, hoisted_automorphism_to!
import ..Rlwe: RLWEkey, SKEncryptor, PKEncryptor, public_keygen, relin_keygen, automorphism_keygen

include("scheme.jl")
include("matrix.jl")
include("polynomial.jl")
include("bootstrap.jl")

export HEScheme, HEOperator, HECiphertext, HEParamSketch, HEParameters, set_encryptor!, encrypt, encrypt_to!, decrypt, decrypt_to!, 
    set_relinkey!, rotate_keygen, set_rotate_key!, set_automorphism_key!, encode, encode_to!, decode, decode_to!, change_level, change_level_to!, drop_level_to!, 
    rescale, rescale_to!, rotate, rotate_to!, hoisted_rotate, hoisted_rotate_to!, get_error
export PlainMatrix, get_bsgs_param, get_required_key_list
export PolyCoeffs, evaluate
export HEBootKey, HEBootstrapper, HEBootParamSketch, HEBootParameters, bootstrap_keygen, set_bootstrapper!, bootstrap

end