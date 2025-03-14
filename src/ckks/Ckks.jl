module Ckks

import ..Math: Modulus, iMform, Bred, Bmul_to!, add, add_to!, sub, sub_to!, neg, neg_to!, gen_power_modmul
import ..Math: UniformSampler, BasisExtender, basis_extend_to!
import ..Math: approximate

import ..Ring: RingParam, find_prime
import ..Ring: PolyEvaluatorNTT, Bred_to!, mul, mul_to!, muladd_to!, ntt!, ntt_to!, intt!, intt_to!,
    automorphism, automorphism_to!, automorphism!, to_big, initialise!
import ..Ring: ComplexPacker, ComplexPackerPow2, pack_to!, unpack_to!
import ..Ring: RefBool, RefInt, RefBigFloat, ComplexBF

import ..Rlwe: PlainConst, PlainPoly, PlainText, RLWE, Tensor, RLEV
import ..Rlwe: decompose, decompose_to!
import ..Rlwe: RLWEParameters
import ..Rlwe: Operator, geteval_at, divide_by_P, divide_by_P_to!, scale, scale_to!, change_modulus, change_modulus_to!, gadgetprod, gadgetprod_to!,
    hoisted_gadgetprod_to!, relinearise, relinearise_to!, keyswitch, keyswitch_to!, hoisted_keyswitch, hoisted_keyswitch_to!,
    hoisted_automorphism, hoisted_automorphism_to!
import ..Rlwe: RLWEkeyPQ, binary_ringkey, ternary_ringkey
import ..Rlwe: SKEncryptor, PKEncryptor, Encryptor, rlwe_sample_to!, rlev_encrypt, automorphism_keygen, relin_keygen, phase, phase_to!

import ..LeveledScheme: HEScheme, HEOperator, HECiphertext, HEParamSketch, HEParameters, encrypt, encrypt_to!, decrypt, decrypt_to!,
    encode, encode_to!, decode, decode_to!, change_level, change_level_to!, drop_level_to!, rescale, rescale_to!, rotate, rotate_to!, hoisted_rotate, hoisted_rotate_to!
import ..LeveledScheme: PlainMatrix, get_bsgs_param, _mul_RLWE!
import ..LeveledScheme: PolyCoeffs, degree, evaluate
import ..LeveledScheme: HEBootKey, bootstrap_keygen, HEBootstrapper, HEBootParamSketch, HEBootParameters, bootstrap

include("ciphertext.jl")
include("parameters.jl")
include("operator.jl")
include("bootstrap.jl")
include("scheme.jl")

export CKKS, CKKSParamSketch, CKKSParameters, CKKSOperator, CKKSScheme
export CKKSBootKey, bootstrap_keygen, CKKSBootParameters, CKKSBootstrapper

end