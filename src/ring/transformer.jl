"""
NTTransformer is an abstract type which supports number theoretic transformations (NTT).
"""
abstract type NTTransformer end

include("./transform/cyclic.jl")
include("./transform/subring.jl")
include("./transform/polybarrett.jl")
include("./transform/cyclotomic.jl")

(::Type{NTTransformer})(param::RingParam, Q::Modulus)::NTTransformer = begin
    if isa(param, CyclicParam)
        CyclicNTTransformer(param.m, Q)     
    elseif isa(param, CyclotomicParam)
        CyclotomicNTTransformer(param.m, Q)
    elseif isa(param, SubringParam)
        SubringNTTransformer(param.m, param.d, Q)
    end
end