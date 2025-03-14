module Ring

import ..Math: Modulus, Moduli, Bred, Bred!, Bred_to!, Bmul, lazy_Bmul, Bmul_to!, lazy_Bmul_to!, Bmuladd_to!, Bmulsub_to!,
    Mform!, iMform!, lazy_Mmul, add, add_to!, sub, sub_to!, neg, neg_to!
import ..Math: BasisExtender, basis_extend_to!
import ..Math: ord, zeropadto, is2a3b5c7d, next2a3b5c7d, factor2357, scramble!, primitive_root_finder,
    ith_root_from_primitive_root, find_generators_mod_m, division, division_mod_Q, cyclotomic_finder
import ..Math: UniformSampler, uniform_random_to!

using Primes: totient, isprime, prevprime, nextprime, factor
using Nemo: finite_field, polynomial_ring, defining_polynomial, ZZ, coeff, lift, zzModPolyRingElem, residue_ring, cyclotomic
using NormalForms: snf
using LinearAlgebra: inv

const ComplexBF = Complex{BigFloat}
const RefBool = Base.RefValue{Bool}
const RefInt = Base.RefValue{UInt64}
const RefBigFloat = Base.RefValue{BigFloat}

include("parameters.jl")
include("transformer.jl")
include("automorphism.jl")
include("modpoly.jl")
include("packing.jl")

export CyclicParam, CyclotomicParam, SubringParam, RingParam, check_modulus, find_prime
export NTTransformer, CyclicNTTransformer, SubringNTTransformer, CyclotomicNTTransformer
export PolyEvaluatorNTT, PolyEvaluatorArb, PolyEvaluatorWord, PolyEvaluatorRNS, ModScalar, ModPoly, Bred, Bred!, Bred_to!,
    mul, mul_to!, muladd_to!, mulsub_to!, reduce!, ntt, ntt!, ntt_to!, intt, intt!, intt_to!, 
    automorphism, automorphism_to!, automorphism!, to_big, uniform_random_to!, initialise!
export ComplexPacker, ComplexPackerPow2, IntPacker, IntPackerArb, IntPackerNTT, IntPackerPow2, IntPackerSubring, pack_to!, unpack_to!
export RefBool, RefInt, RefBigFloat, ComplexBF

end