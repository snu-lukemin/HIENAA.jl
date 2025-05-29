function _automorphism!(idx::Int64, a::AbstractVector{UInt64}, isntt::Bool, ntter::SubringNTTransformer)::Nothing
    m, gpowN, ginvpowN = ntter.m, ntter.gpowN, ntter.ginvpowN

    idx = mod(idx, m)

    if idx ∈ gpowN
        shift = findfirst(isequal(idx), gpowN) - 1
    elseif idx ∈ ginvpowN
        shift = -(findfirst(isequal(idx), ginvpowN) - 1)
    else
        throw(DomainError("$(idx) is not a valid automorphism index."))
    end

    circshift!(a, isntt ? -shift : shift)

    return nothing
end