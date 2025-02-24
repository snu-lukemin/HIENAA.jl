"""
    HIENAA

A package for RLWE-based homomorphic encryption schemes.
"""
module HIENAA

using ChaChaCiphers: ChaChaStream, ChaCha12Stream
using Primes: totient, isprime, prevprime, nextprime, factor
using Nemo: finite_field, polynomial_ring, defining_polynomial, ZZ, coeff, lift, zzModPolyRingElem, residue_ring, cyclotomic
using NormalForms: snf
using LinearAlgebra: inv
using DoubleFloats: Double64, ComplexDF64
export Double64, ComplexDF64

const RefBool = Base.RefValue{Bool}
const RefInt = Base.RefValue{UInt64}
const RefDouble = Base.RefValue{Double64}

include("math/modular.jl")
include("math/basischange.jl")
include("math/arithmetic.jl")
include("math/sampler.jl")
include("math/polynomials.jl")

include("ring/parameters.jl")
include("ring/transformer.jl")
include("ring/automorphism.jl")
include("ring/modpoly.jl")
include("ring/packing.jl")

include("rlwe/ciphertext.jl")
include("rlwe/decomposition.jl")
include("rlwe/parameters.jl")
include("rlwe/operator.jl")
include("rlwe/key.jl")
include("rlwe/encryptor.jl")

include("bgv/ciphertext.jl")
include("bgv/parameters.jl")
include("bgv/operator.jl")

include("bfv/ciphertext.jl")
include("bfv/parameters.jl")
include("bfv/operator.jl")

include("ckks/ciphertext.jl")
include("ckks/parameters.jl")
include("ckks/operator.jl")

include("common/scheme.jl")
include("common/polynomial.jl")
include("common/matrix.jl")

export encode, encode_to!, decode, decode_to!, change_level, change_level_to!, drop_level, drop_level_to!,
       rescale, rescale_to!, neg, neg_to!, add, add_to!, sub, sub_to!, mul, mul_to!, mul, mul_to!, mul, mul_to!,
       keyswitch, keyswitch_to!, hoisted_keyswitch, hoisted_keyswitch_to!, rotate, rotate_to!, hoisted_rotate, hoisted_rotate_to!, 
       automorphism, automorphism_to!, hoisted_automorphism, hoisted_automorphism_to!, decompose_a, decompose_a_to!

end