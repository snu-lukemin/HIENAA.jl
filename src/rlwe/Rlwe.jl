module Rlwe

import ..Math: Modulus, Moduli, Bred, Bmul, add, add_to!, sub, sub_to!, neg, neg_to!
import ..Math: BasisExtender, basis_extend_to!, SimpleScaler, simple_scale_to!
import ..Math: UniformSampler, uniform_binary, uniform_ternary, uniform_random, uniform_random_to!,
    CDTSampler, VarCenSampler, TwinCDTSampler, RGSampler, sample

import ..Ring: CyclicParam, CyclotomicParam, SubringParam, RingParam, check_modulus, find_prime
import ..Ring: PolyEvaluatorNTT, PolyEvaluatorArb, PolyEvaluatorWord, PolyEvaluatorRNS, ModScalar, ModPoly, mul, mul_to!,
    muladd_to!, mulsub_to!, reduce!, ntt, ntt!, ntt_to!, intt, intt!, intt_to!, automorphism, automorphism_to!, automorphism!, 
    to_big, initialise!
import ..Ring: RefBool, RefInt, RefBigFloat

include("ciphertext.jl")
include("decomposition.jl")
include("parameters.jl")
include("operator.jl")
include("key.jl")
include("encryptor.jl")
include("extraction.jl")

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

end