module Bfv

import ..Math: Modulus, Bred, Bmul, Bmul_to!, add, add_to!, sub, sub_to!, neg, neg_to!, gen_power_modmul
import ..Math: BasisExtender, basis_extend_to!, SimpleScaler, simple_scale_to!, ComplexScaler, complex_scale_to!
import ..Math: Interpolator, interpolate

import ..Ring: RingParam, find_prime
import ..Ring: PolyEvaluatorRNS, ModPoly, Bred_to!, mul, mul_to!, muladd_to!, mulsub_to!, ntt, ntt!, ntt_to!, intt, intt!, intt_to!,
    automorphism, automorphism_to!, automorphism!, to_big
import ..Ring: IntPacker, IntPackerArb, IntPackerNTT, IntPackerPow2, IntPackerSubring, pack_to!, unpack_to!, load_resolution
import ..Ring: RefBool, RefInt

import ..Rlwe: PlainConst, PlainPoly, PlainText, RLWE, Tensor, RLEV, RGSW
import ..Rlwe: Decomposer, decompose, decompose_to!
import ..Rlwe: RLWEParamSketch, RLWEParameters
import ..Rlwe: Operator, geteval_at, scale, scale_to!, change_modulus_to!, relinearise_to!, keyswitch, keyswitch_to!,
    hoisted_keyswitch, hoisted_keyswitch_to!, hoisted_automorphism, hoisted_automorphism_to!
import ..Rlwe: SKEncryptor, PKEncryptor, Encryptor, rlwe_sample_to!, phase, phase_to!
import ..Rlwe: Extractor, ExtractorSubring, extract

import ..LeveledScheme: HEScheme, HEOperator, HECiphertext, HEParamSketch, HEParameters, set_encryptor!, encrypt, encrypt_to!, decrypt, decrypt_to!,
    set_relinkey!, rotate_keygen, set_rotate_key!, set_automorphism_key!, encode, encode_to!, decode, decode_to!, change_level, change_level_to!, drop_level_to!,
    rescale, rescale_to!, rotate, rotate_to!, hoisted_rotate, hoisted_rotate_to!, get_error
import ..LeveledScheme: PlainMatrix, get_bsgs_param, _mul_RLWE!
import ..LeveledScheme: PolyCoeffs, degree, evaluate
import ..LeveledScheme: HEBootstrapper, HEBootParamSketch, HEBootParameters, bootstrap

include("ciphertext.jl")
include("parameters.jl")
include("operator.jl")
include("scheme.jl")
include("bootstrap.jl")

export BFV, BFVParamSketch, BFVParameters, BFVOperator, BFVScheme

end