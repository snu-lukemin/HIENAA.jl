function _automorphism!(idx::Int64, a::AbstractVector{UInt64}, isntt::Bool, ntter::SubringNTTransformer)
    m, d, N, g, gpowN = ntter.m, ntter.d, ntter.N, ntter.g, ntter.gpowN

    idx = mod(idx, m)

    gᴺ = powermod(g, N, m)
    for i = 0:d-1
        tmp = (idx * powermod(gᴺ, i, m)) % m
        if tmp ∈ gpowN
            shift = findfirst(isequal(tmp), gpowN) - 1
            circshift!(a, isntt ? -shift : shift)
        end
    end
end