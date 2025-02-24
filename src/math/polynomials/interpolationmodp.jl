function ν(m::Int64, p::Int64)
    res = 0
    while m % p == 0
        m ÷= p
        res += 1
    end
    res
end

function μ(p::Int64, e::Int64)
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

    function Interpolator(p::Int64, r::Int64)
        pr = Modulus(p^r)
        buff = [zeros(UInt64, μ(p, r)) for _ = 1:5]
        new(p, r, pr, buff)
    end
end

@views function finite_diff_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, interp::Interpolator)
    len = length(x)
    pr, p, r = interp.pr, interp.p, interp.r

    @assert len ≤ μ(p, r) "The number of points is too large."
    @assert length(res) == length(y) == len "The length of the input vectors must be the same."

    tmp1, tmp2 = interp.buff[4][1:len], interp.buff[5][1:len]
    
    @. res = 0
    @. tmp1 = y
    @. tmp2 = 0
    @inbounds for i = 2:len
        res[i-1] = tmp1[1]
        for j = 1:len-i+1
            resj = _sub(tmp1[j+1], tmp1[j], pr)
            denom = _sub(x[j+i-1], x[j], pr)

            while denom % p == 0
                denom ÷= p
                @assert resj % p == 0 "The input function is not a polynomial."
                resj ÷= p
            end

            tmp2[j] = _Bmul(resj, invmod(denom, pr.Q), pr)
        end
        @. tmp1 = tmp2
    end
    res[end] = tmp1[1]

    res
end

@views function interpolate(xi::AbstractVector{<:Integer}, yi::AbstractVector{<:Integer}, interp::Interpolator)
    pr, p, r = interp.pr, interp.p, interp.r

    @assert length(xi) == length(yi) "The length of the input vectors must be the same."
    @assert length(xi) ≤ μ(p, r) "The number of points is too large."

    len = length(xi)
    x, y, diff = interp.buff[1][1:len], interp.buff[2][1:len], interp.buff[3][1:len]
    
    _Bred_to!(x, xi, pr)
    _Bred_to!(y, yi, pr)

    finite_diff_to!(diff, x, y, interp)
    res = zeros(UInt64, len)

    basis = interp.buff[4][1:len]
    @. basis = 0
    basis[1] = 1
    @inbounds for i = 1:len
        _Bmuladd_to!(res, diff[i], basis, pr)

        if i < len
            @. basis[2:i+1] = basis[1:i]
            basis[1] = 0
            _Bmuladd_to!(basis[1:i], pr.Q - x[i], basis[2:i+1], pr)
        end
    end

    res
end

export Interpolator, interpolate