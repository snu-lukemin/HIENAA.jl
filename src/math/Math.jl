"""
    HIENAA.Math

An internal submodule for math.
"""
module Math

using ChaChaCiphers: ChaChaStream, ChaCha12Stream
using Primes:totient, factor

include("modular.jl")
include("basischange.jl")
include("arithmetic.jl")
include("sampler.jl")
include("polynomials.jl")

export Modulus, Moduli, Bred, Bred, lazy_Bred, Bred!, Bred_to!, Bmul, lazy_Bmul, Bmul_to!, lazy_Bmul_to!, 
    Mform, iMform, Mform!, Mform_to!, iMform!, iMform_to!, Mmul, lazy_Mmul, Mmul_to!, 
    add, add_to!, sub, sub_to!, neg, neg_to!, Bmuladd_to!, Bmulsub_to!
export BasisExtender, basis_extend_to!, SimpleScaler, simple_scale_to!, ComplexScaler, complex_scale_to!
export ord, zeropadto, is2a3b5c7d, next2a3b, factor2357, mult_and_round, round_to_uint64, round_to_uint128,
    scramble!, primitive_root_finder, ith_root_from_primitive_root, find_generators_mod_m, division, division_mod_Q,
    cyclotomic_finder, gen_power_modmul, round
export UniformSampler, uniform_binary, uniform_ternary, uniform_random, uniform_random_to!,
    CDTSampler, VarCenSampler, TwinCDTSampler, RGSampler, sample
export Interpolator, interpolate
export approximate

end