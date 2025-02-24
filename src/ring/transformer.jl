include("./transform/cyclic.jl")
include("./transform/subring.jl")
include("./transform/polybarrett.jl")
include("./transform/cyclotomic.jl")

"""
NTTransformer is an abstract type which supports number theoretic transformations (NTT).
"""
const NTTransformer = Union{CyclicNTTransformer,SubringNTTransformer,CyclotomicNTTransformer}

(::Type{NTTransformer})(param::RingParam, Q::Modulus) = begin
    if typeof(param) == CyclicParam
        CyclicNTTransformer(param.m, Q)     
    elseif typeof(param) == CyclotomicParam
        CyclotomicNTTransformer(param.m, Q)
    elseif typeof(param) == SubringParam
        SubringNTTransformer(param.m, param.d, Q)
    end
end

export CyclicNTTransformer, SubringNTTransformer, CyclotomicNTTransformer, NTTransformer