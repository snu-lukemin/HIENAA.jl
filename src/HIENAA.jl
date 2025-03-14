"""
    HIENAA

A package for RLWE-based homomorphic encryption schemes.
"""
module HIENAA

include("math/Math.jl")
include("ring/Ring.jl")
include("rlwe/Rlwe.jl")
include("common/LeveledScheme.jl")
include("bgv/Bgv.jl")
include("bfv/Bfv.jl")
include("ckks/Ckks.jl")

using .Math
using .Ring
using .Rlwe
using .LeveledScheme
using .Bgv
using .Bfv
using .Ckks

# Double
export ComplexBF

# Math
export Modulus, Moduli, BasisExtender, basis_extend_to!, SimpleScaler, simple_scale_to!, ComplexScaler, complex_scale_to!,
    UniformSampler, uniform_binary, uniform_ternary, uniform_random, uniform_random_to!,
    CDTSampler, VarCenSampler, TwinCDTSampler, RGSampler, sample, Interpolator, interpolate

# Ring
export CyclicParam, CyclotomicParam, SubringParam, RingParam, check_modulus, find_prime,
    NTTransformer, CyclicNTTransformer, SubringNTTransformer, CyclotomicNTTransformer,
    PolyEvaluatorNTT, PolyEvaluatorArb, PolyEvaluator, PolyEvaluatorRNS, ModScalar, ModPoly, Bred, Bred!, Bred_to!,
    neg, neg_to!, add, add_to!, sub, sub_to!, mul, mul_to!, muladd_to!, mulsub_to!, reduce!, ntt, ntt!, ntt_to!, intt, intt!, intt_to!,
    automorphism, automorphism_to!, automorphism!, to_big, uniform_random_to!, ComplexPacker, IntPacker, pack_to!, unpack_to!,
    RefBool, RefInt, RefBigFloat

# Rlwe
export PlainConst, PlainPoly, PlainText, RLWE, Tensor, RLEV, RGSW
export Decomposer, decompose, decompose_to!
export RLWEParamSketch, RLWEParameters
export Operator, geteval_at, divide_by_P, divide_by_P_to!, scale, scale_to!, change_modulus, change_modulus_to!, gadgetprod, gadgetprod_to!,
    hoisted_gadgetprod_to!, relinearise, relinearise_to!, keyswitch, keyswitch_to!, hoisted_keyswitch, hoisted_keyswitch_to!,
    hoisted_automorphism, hoisted_automorphism_to!, extprod, extprod_to!, hoisted_extprod_to!
export RLWEkey, RLWEkeyPQ, binary_ringkey, ternary_ringkey
export SKEncryptor, PKEncryptor, Encryptor, rlwe_sample, rlwe_sample_to!, rlwe_encrypt, rlwe_encrypt_to!, rlwe_encrypt_a, rlwe_encrypt_a_to!,
    phase, phase_to!, rlev_encrypt, rlev_encrypt_to!, rlev_encrypt_a, rlev_encrypt_a_to!, rgsw_encrypt, rgsw_encrypt_to!,
    public_keygen, relin_keygen, automorphism_keygen
export Extractor, ExtractorSubring, extract

# LeveledScheme
export HEScheme, HEOperator, HECiphertext, HEParamSketch, HEParameters, set_encryptor!, encrypt, encrypt_to!, decrypt, decrypt_to!,
    set_relinkey!, rotate_keygen, set_rotate_key!, set_automorphism_key!, encode, encode_to!, decode, decode_to!, change_level, change_level_to!, 
    drop_level, drop_level_to!, rescale, rescale_to!, rotate, rotate_to!, hoisted_rotate, hoisted_rotate_to!
export PlainMatrix, get_bsgs_param, get_required_key_list
export PolyCoeffs, degree, evaluate
export HEBootstrapper, HEBootParamSketch, HEBootParameters, bootstrap_keygen, set_bootstrapper!, bootstrap

# Bgv
export BGV, BGVParamSketch, BGVParameters, BGVOperator, BGVScheme

# Bfv
export BFV, BFVParamSketch, BFVParameters, BFVOperator, BFVScheme

# Ckks
export CKKS, CKKSParamSketch, CKKSParameters, CKKSOperator, CKKSScheme, CKKSBootParameters, CKKSBootKey, CKKSBootstrapper

end