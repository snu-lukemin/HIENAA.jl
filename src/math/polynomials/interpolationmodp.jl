function ν(m::Int64, p::Int64)::Int64
    res = 0
    while m % p == 0
        m ÷= p
        res += 1
    end
    res
end

function μ(p::Int64, e::Int64)::Int64
    i, cnt = p, 0
    while cnt < e
        cnt += ν(i, p)
        i += p
    end
    i - p
end

struct Interpolator
    p::Int64
    r::Int64
    pr::Modulus
    buff::Vector{Vector{UInt64}}

    function Interpolator(p::Int64, r::Int64)::Interpolator
        pr = Modulus(p^r)
        buff = [zeros(UInt64, μ(p, r)) for _ = 1:5]
        new(p, r, pr, buff)
    end
end

@views function finite_diff_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, interp::Interpolator)::Nothing
    len = length(x)
    pr, p, r = interp.pr, interp.p, interp.r

    if len > μ(p, r)
        throw(DomainError("The number of points is too large."))
    end
    if length(res) ≠ length(y) ≠ len
        throw(DimensionMismatch("The length of the input vectors must be the same."))
    end

    tmp1, tmp2 = interp.buff[4][1:len], interp.buff[5][1:len]
    
    @. res = 0
    @. tmp1 = y
    @. tmp2 = 0
    @inbounds for i = 2:len
        res[i-1] = tmp1[1]
        for j = 1:len-i+1
            resj = sub(tmp1[j+1], tmp1[j], pr)
            denom = sub(x[j+i-1], x[j], pr)

            while denom % p == 0
                denom ÷= p
                if resj % p ≠ 0
                    throw(DomainError("The input function is not a polynomial."))
                end
                resj ÷= p
            end

            tmp2[j] = Bmul(resj, invmod(denom, pr.Q), pr)
        end
        @. tmp1 = tmp2
    end
    res[end] = tmp1[1]

    return nothing
end

@views function interpolate(xi::AbstractVector{<:Integer}, yi::AbstractVector{<:Integer}, interp::Interpolator)::Vector{UInt64}
    pr, p, r = interp.pr, interp.p, interp.r

    if length(xi) ≠ length(yi)
        throw(DimensionMismatch("The length of the input vectors must be the same."))
    end
    if length(xi) > μ(p, r)
        throw(DomainError("The number of points is too large."))
    end

    len = length(xi)
    x, y, diff = interp.buff[1][1:len], interp.buff[2][1:len], interp.buff[3][1:len]
    
    Bred_to!(x, xi, pr)
    Bred_to!(y, yi, pr)

    finite_diff_to!(diff, x, y, interp)
    res = zeros(UInt64, len)

    basis = interp.buff[4][1:len]
    @. basis = 0
    basis[1] = 1
    @inbounds for i = 1:len
        Bmuladd_to!(res, diff[i], basis, pr)

        if i < len
            @. basis[2:i+1] = basis[1:i]
            basis[1] = 0
            Bmuladd_to!(basis[1:i], pr.Q - x[i], basis[2:i+1], pr)
        end
    end

    res
end