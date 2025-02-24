"""
RLWEkey is a struct for the RLWE secret key s ∈ R.
"""
struct RLWEkey
    coeffs::Vector{Int64}
    N::Int64
end

"""
Outputs secret key for RLWE.
If hw = 0, the secret is sampled from a uniform binary distribution.
Otherwise, it outputs a secret key with hamming weight at most hw. 
"""
binary_ringkey(us::UniformSampler, N::Int64, hw::Int64=0) = RLWEkey(uniform_binary(us, N, hw), N)

"""
Outputs secret key for RLWE.
If hw = 0, the secret is sampled from a uniform ternary distribution.
Otherwise, it outputs a secret key with hamming weight at most hw. 
"""
ternary_ringkey(us::UniformSampler, N::Int64, hw::Int64=0) = RLWEkey(uniform_ternary(us, N, hw), N)

"""
RLWEkeyPQ is a struct for the embedding of the RLWE secret key s ∈ R_Q.
"""
const RLWEkeyPQ = ModPoly

RLWEkeyPQ(key::RLWEkey, eval::PolyEvaluator) = begin
    coeffs = [Vector{UInt64}(undef, key.N) for _ = eachindex(eval)]
    @inbounds for i = eachindex(eval)
        for j = 1:key.N
            coeffs[i][j] = _Bred(key.coeffs[j], eval[i].Q)
        end
    end
    ModPoly(coeffs, false)
end

export RLWEkey, RLWEkeyPQ, binary_ringkey, ternary_ringkey