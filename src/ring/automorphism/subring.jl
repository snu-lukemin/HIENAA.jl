function _automorphism!(idx::Int64, a::AbstractVector{UInt64}, isntt::Bool, ntter::SubringNTTransformer)::Nothing
    m, d, N, g, gpowN = ntter.m, ntter.d, ntter.N, ntter.g, ntter.gpowN

    idx = mod(idx, m)

    if idx ∉ gpowN
        throw(DomainError("$(idx) is not a valid automorphism index."))
    end

    shift = findfirst(isequal(idx), gpowN) - 1
    circshift!(a, isntt ? shift : -shift)

    return nothing
end