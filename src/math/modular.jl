"""
Modulus supports a modular arithmetic given modulus Q.
"""
struct Modulus
    Q::UInt64
    Q⁻¹::UInt64       # Used for Montgomery reduction.
    R⁻¹::UInt64       # Used for Montgomery reduction.
    r0::UInt64        # Used for Barrett reduction.
    r1::UInt64        # Used for Barrett reduction.
    logQ::Int64       # Used for gadget decomposition.
    halfQ::Int64      # Used for balanced representation.

    function Modulus(Q::Integer)
        Q == 0 && return new(0, 0, 0, 0, 0, 0, 0)

        logQ = ceil(Int64, log2(Q))
        @assert logQ ≤ 62 "Modulus should be smaller than 2⁶²."

        if isodd(Q)
            Q⁻¹ = UInt64(invmod(-Q, 0x00000000000000010000000000000000))
            R⁻¹ = UInt64(invmod(0x00000000000000010000000000000000, Q))
        else
            Q⁻¹ = UInt64(0)
            R⁻¹ = UInt128(0)
        end

        r = floor(UInt128, (big(1) << 128) / Q)
        r0, r1 = r % UInt64, (r >> 64) % UInt64
        halfQ = (Q - 1) ÷ 2
        new(Q, Q⁻¹, R⁻¹, r0, r1, logQ, halfQ)
    end
end

Base.:show(io::IO, Q::Modulus) = print(io, Q.Q)
Base.:show(io::IO, mime::MIME"text/plain", Q::Modulus) = println(io, Q.Q)

"""
_Bred(x, Q) returns x % Q using barret reduction.
"""
_Bred(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) >> 64) % UInt64
    res = x - Q.Q * t
    res ≥ Q.Q ? res - Q.Q : res
end

_Bred(x::Int64, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return _Bred(UInt64(x), Q)
    else
        res = _Bred(UInt64(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

_Bred(x::UInt128, Q::Modulus)::UInt64 = begin
    x0 = x % UInt64
    x1 = (x >> 64) % UInt64
    t = (widemul(Q.r1, x1) + (widemul(Q.r1, x0) + widemul(Q.r0, x1) + widemul(Q.r0, x0) >> 64) >> 64) % UInt64
    res = x0 - t * Q.Q
    res ≥ Q.Q ? res - Q.Q : res
end

_Bred(x::Int128, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return _Bred(UInt128(x), Q)
    else
        res = _Bred(UInt128(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

_lazy_Bred(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) >> 64) % UInt64
    x - Q.Q * t
end

_lazy_Bred(x::Int64, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return _lazy_Bred(UInt64(x), Q)
    else
        res = _lazy_Bred(UInt64(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

_lazy_Bred(x::UInt128, Q::Modulus)::UInt64 = begin
    x0 = x % UInt64
    x1 = (x >> 64) % UInt64
    t = (widemul(Q.r1, x1) + (widemul(Q.r1, x0) + widemul(Q.r0, x1) + widemul(Q.r0, x0) >> 64) >> 64) % UInt64
    x0 - t * Q.Q
end

_lazy_Bred(x::Int128, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return _lazy_Bred(UInt128(x), Q)
    else
        res = _lazy_Bred(UInt128(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

_Bred!(a::AbstractVector{UInt64}, Q::Modulus) = _Bred_to!(a, a, Q)

_Bred_to!(res::AbstractVector{UInt64}, a::AbstractVector{<:Union{Int64,UInt64,Int128,UInt128}}, Q::Modulus) = begin
    @assert length(res) == length(a) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bred(a[i+1], Q)
        res[i+2] = _Bred(a[i+2], Q)
        res[i+3] = _Bred(a[i+3], Q)
        res[i+4] = _Bred(a[i+4], Q)
        res[i+5] = _Bred(a[i+5], Q)
        res[i+6] = _Bred(a[i+6], Q)
        res[i+7] = _Bred(a[i+7], Q)
        res[i+8] = _Bred(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bred(a[i], Q)
    end
end

# Taken from the invmod.
Base.:invmod(x::Integer, Q::Modulus)::UInt64 = begin
    n, m = UInt64(x), Q.Q
    g, x, _ = gcdx(n, m)
    g ≠ 1 && throw(DomainError((n, m), LazyString("Greatest common divisor is ", g, ".")))
    x > typemax(UInt64) >> 1 && (x += m)
    _Bred(x, Q)
end

_Mform(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r0, x) >> 64) % UInt64 + Q.r1 * x
    res = -t * Q.Q
    res ≥ Q.Q ? res - Q.Q : res
end

_Mform(x::Int64, Q::Modulus)::UInt64 = _Mform(Bred(x, Q), Q)

_iMform(x::UInt64, Q::Modulus)::UInt64 = _Bred(widemul(x, Q.R⁻¹), Q)

_Mform!(x::AbstractVector{UInt64}, Q::Modulus) = _Mform_to!(x, x, Q)

_Mform_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Mform(a[i+1], Q)
        res[i+2] = _Mform(a[i+2], Q)
        res[i+3] = _Mform(a[i+3], Q)
        res[i+4] = _Mform(a[i+4], Q)
        res[i+5] = _Mform(a[i+5], Q)
        res[i+6] = _Mform(a[i+6], Q)
        res[i+7] = _Mform(a[i+7], Q)
        res[i+8] = _Mform(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Mform(a[i], Q)
    end
end

_iMform!(x::AbstractVector{UInt64}, Q::Modulus) = _iMform_to!(x, x, Q)

_iMform_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _iMform(a[i+1], Q)
        res[i+2] = _iMform(a[i+2], Q)
        res[i+3] = _iMform(a[i+3], Q)
        res[i+4] = _iMform(a[i+4], Q)
        res[i+5] = _iMform(a[i+5], Q)
        res[i+6] = _iMform(a[i+6], Q)
        res[i+7] = _iMform(a[i+7], Q)
        res[i+8] = _iMform(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _iMform(a[i], Q)
    end
end

"""
Convert UInt64 to a balanced representation.
"""
Base.:signed(x::UInt64, Q::Modulus)::Int64 = signed(x > Q.halfQ ? x - Q.Q : x)

_add(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    res = a + b
    res ≥ Q.Q ? res - Q.Q : res
end

_add_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _add(a[i+1], b[i+1], Q)
        res[i+2] = _add(a[i+2], b[i+2], Q)
        res[i+3] = _add(a[i+3], b[i+3], Q)
        res[i+4] = _add(a[i+4], b[i+4], Q)
        res[i+5] = _add(a[i+5], b[i+5], Q)
        res[i+6] = _add(a[i+6], b[i+6], Q)
        res[i+7] = _add(a[i+7], b[i+7], Q)
        res[i+8] = _add(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _add(a[i], b[i], Q)
    end
end

_add_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::UInt64, Q::Modulus) = begin
    @assert length(res) == length(a) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _add(a[i+1], b, Q)
        res[i+2] = _add(a[i+2], b, Q)
        res[i+3] = _add(a[i+3], b, Q)
        res[i+4] = _add(a[i+4], b, Q)
        res[i+5] = _add(a[i+5], b, Q)
        res[i+6] = _add(a[i+6], b, Q)
        res[i+7] = _add(a[i+7], b, Q)
        res[i+8] = _add(a[i+8], b, Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _add(a[i], b, Q)
    end
end

_neg(a::UInt64, Q::Modulus)::UInt64 = a == 0 ? a : Q.Q - a

_neg_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _neg(a[i+1], Q)
        res[i+2] = _neg(a[i+2], Q)
        res[i+3] = _neg(a[i+3], Q)
        res[i+4] = _neg(a[i+4], Q)
        res[i+5] = _neg(a[i+5], Q)
        res[i+6] = _neg(a[i+6], Q)
        res[i+7] = _neg(a[i+7], Q)
        res[i+8] = _neg(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _neg(a[i], Q)
    end
end

_sub(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    res = a - b
    res ≥ Q.Q ? res + Q.Q : res
end

_sub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _sub(a[i+1], b[i+1], Q)
        res[i+2] = _sub(a[i+2], b[i+2], Q)
        res[i+3] = _sub(a[i+3], b[i+3], Q)
        res[i+4] = _sub(a[i+4], b[i+4], Q)
        res[i+5] = _sub(a[i+5], b[i+5], Q)
        res[i+6] = _sub(a[i+6], b[i+6], Q)
        res[i+7] = _sub(a[i+7], b[i+7], Q)
        res[i+8] = _sub(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _sub(a[i], b[i], Q)
    end
end

_sub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::UInt64, Q::Modulus) = begin
    @assert length(res) == length(a) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _sub(a[i+1], b, Q)
        res[i+2] = _sub(a[i+2], b, Q)
        res[i+3] = _sub(a[i+3], b, Q)
        res[i+4] = _sub(a[i+4], b, Q)
        res[i+5] = _sub(a[i+5], b, Q)
        res[i+6] = _sub(a[i+6], b, Q)
        res[i+7] = _sub(a[i+7], b, Q)
        res[i+8] = _sub(a[i+8], b, Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _sub(a[i], b, Q)
    end
end

_Mmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    q, q⁻¹ = Q.Q, Q.Q⁻¹

    ab = widemul(a, b)
    w = ((widemul(q, (ab % UInt64) * q⁻¹) + ab) >> 64) % UInt64
    w ≥ q ? w - q : w
end

_lazy_Mmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    q, q⁻¹ = Q.Q, Q.Q⁻¹

    ab = widemul(a, b)
    ((widemul(q, (ab % UInt64) * q⁻¹) + ab) >> 64) % UInt64
end

_Bmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = _Bred(widemul(a, b), Q)

_lazy_Bmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = _lazy_Bred(widemul(a, b), Q)

_Mmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Mmul(a[i+1], b[i+1], Q)
        res[i+2] = _Mmul(a[i+2], b[i+2], Q)
        res[i+3] = _Mmul(a[i+3], b[i+3], Q)
        res[i+4] = _Mmul(a[i+4], b[i+4], Q)
        res[i+5] = _Mmul(a[i+5], b[i+5], Q)
        res[i+6] = _Mmul(a[i+6], b[i+6], Q)
        res[i+7] = _Mmul(a[i+7], b[i+7], Q)
        res[i+8] = _Mmul(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Mmul(a[i], b[i], Q)
    end
end

_Mmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(b) "Length of res and a should be the same."
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Mmul(a, b[i+1], Q)
        res[i+2] = _Mmul(a, b[i+2], Q)
        res[i+3] = _Mmul(a, b[i+3], Q)
        res[i+4] = _Mmul(a, b[i+4], Q)
        res[i+5] = _Mmul(a, b[i+5], Q)
        res[i+6] = _Mmul(a, b[i+6], Q)
        res[i+7] = _Mmul(a, b[i+7], Q)
        res[i+8] = _Mmul(a, b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Mmul(a, b[i], Q)
    end
end

_Bmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bmul(a[i+1], b[i+1], Q)
        res[i+2] = _Bmul(a[i+2], b[i+2], Q)
        res[i+3] = _Bmul(a[i+3], b[i+3], Q)
        res[i+4] = _Bmul(a[i+4], b[i+4], Q)
        res[i+5] = _Bmul(a[i+5], b[i+5], Q)
        res[i+6] = _Bmul(a[i+6], b[i+6], Q)
        res[i+7] = _Bmul(a[i+7], b[i+7], Q)
        res[i+8] = _Bmul(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bmul(a[i], b[i], Q)
    end
end

_Bmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(b) "Length of res and a should be the same."
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bmul(a, b[i+1], Q)
        res[i+2] = _Bmul(a, b[i+2], Q)
        res[i+3] = _Bmul(a, b[i+3], Q)
        res[i+4] = _Bmul(a, b[i+4], Q)
        res[i+5] = _Bmul(a, b[i+5], Q)
        res[i+6] = _Bmul(a, b[i+6], Q)
        res[i+7] = _Bmul(a, b[i+7], Q)
        res[i+8] = _Bmul(a, b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bmul(a, b[i], Q)
    end
end

_lazy_Bmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _lazy_Bmul(a[i+1], b[i+1], Q)
        res[i+2] = _lazy_Bmul(a[i+2], b[i+2], Q)
        res[i+3] = _lazy_Bmul(a[i+3], b[i+3], Q)
        res[i+4] = _lazy_Bmul(a[i+4], b[i+4], Q)
        res[i+5] = _lazy_Bmul(a[i+5], b[i+5], Q)
        res[i+6] = _lazy_Bmul(a[i+6], b[i+6], Q)
        res[i+7] = _lazy_Bmul(a[i+7], b[i+7], Q)
        res[i+8] = _lazy_Bmul(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _lazy_Bmul(a[i], b[i], Q)
    end
end

_lazy_Bmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(b) "Length of res and a should be the same."
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = _lazy_Bmul(a, b[i+1], Q)
        res[i+2] = _lazy_Bmul(a, b[i+2], Q)
        res[i+3] = _lazy_Bmul(a, b[i+3], Q)
        res[i+4] = _lazy_Bmul(a, b[i+4], Q)
        res[i+5] = _lazy_Bmul(a, b[i+5], Q)
        res[i+6] = _lazy_Bmul(a, b[i+6], Q)
        res[i+7] = _lazy_Bmul(a, b[i+7], Q)
        res[i+8] = _lazy_Bmul(a, b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _lazy_Bmul(a, b[i], Q)
    end
end

_Bmuladd_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bred(widemul(a[i+1], b[i+1]) + res[i+1], Q)
        res[i+2] = _Bred(widemul(a[i+2], b[i+2]) + res[i+2], Q)
        res[i+3] = _Bred(widemul(a[i+3], b[i+3]) + res[i+3], Q)
        res[i+4] = _Bred(widemul(a[i+4], b[i+4]) + res[i+4], Q)
        res[i+5] = _Bred(widemul(a[i+5], b[i+5]) + res[i+5], Q)
        res[i+6] = _Bred(widemul(a[i+6], b[i+6]) + res[i+6], Q)
        res[i+7] = _Bred(widemul(a[i+7], b[i+7]) + res[i+7], Q)
        res[i+8] = _Bred(widemul(a[i+8], b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bred(widemul(a[i], b[i]) + res[i], Q)
    end
end

_Bmuladd_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(b) "Length of res and a should be the same."
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bred(widemul(a, b[i+1]) + res[i+1], Q)
        res[i+2] = _Bred(widemul(a, b[i+2]) + res[i+2], Q)
        res[i+3] = _Bred(widemul(a, b[i+3]) + res[i+3], Q)
        res[i+4] = _Bred(widemul(a, b[i+4]) + res[i+4], Q)
        res[i+5] = _Bred(widemul(a, b[i+5]) + res[i+5], Q)
        res[i+6] = _Bred(widemul(a, b[i+6]) + res[i+6], Q)
        res[i+7] = _Bred(widemul(a, b[i+7]) + res[i+7], Q)
        res[i+8] = _Bred(widemul(a, b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bred(widemul(a, b[i]) + res[i], Q)
    end
end

_Bmulsub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(a) == length(b) "Length of res and a should be the same."
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bred(widemul(Q.Q - a[i+1], b[i+1]) + res[i+1], Q)
        res[i+2] = _Bred(widemul(Q.Q - a[i+2], b[i+2]) + res[i+2], Q)
        res[i+3] = _Bred(widemul(Q.Q - a[i+3], b[i+3]) + res[i+3], Q)
        res[i+4] = _Bred(widemul(Q.Q - a[i+4], b[i+4]) + res[i+4], Q)
        res[i+5] = _Bred(widemul(Q.Q - a[i+5], b[i+5]) + res[i+5], Q)
        res[i+6] = _Bred(widemul(Q.Q - a[i+6], b[i+6]) + res[i+6], Q)
        res[i+7] = _Bred(widemul(Q.Q - a[i+7], b[i+7]) + res[i+7], Q)
        res[i+8] = _Bred(widemul(Q.Q - a[i+8], b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bred(widemul(Q.Q - a[i], b[i]) + res[i], Q)
    end
end

_Bmulsub_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus) = begin
    @assert length(res) == length(b) "Length of res and a should be the same."
    N = length(b)
    neg_a = Q.Q - a
    @inbounds for i = 0:8:N-8
        res[i+1] = _Bred(widemul(neg_a, b[i+1]) + res[i+1], Q)
        res[i+2] = _Bred(widemul(neg_a, b[i+2]) + res[i+2], Q)
        res[i+3] = _Bred(widemul(neg_a, b[i+3]) + res[i+3], Q)
        res[i+4] = _Bred(widemul(neg_a, b[i+4]) + res[i+4], Q)
        res[i+5] = _Bred(widemul(neg_a, b[i+5]) + res[i+5], Q)
        res[i+6] = _Bred(widemul(neg_a, b[i+6]) + res[i+6], Q)
        res[i+7] = _Bred(widemul(neg_a, b[i+7]) + res[i+7], Q)
        res[i+8] = _Bred(widemul(neg_a, b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = _Bred(widemul(neg_a, b[i]) + res[i], Q)
    end
end

const Moduli = AbstractVector{Modulus}

Base.:prod(Q::Moduli) = begin
    Qbig = BigInt(1)
    @inbounds for i = eachindex(Q)
        Qbig *= Q[i].Q
    end
    Qbig
end

Base.:log2(Q::Moduli) = begin
    logQ = 0.0
    @inbounds for i = eachindex(Q)
        logQ += log2(Q[i].Q)
    end
    logQ
end

export Modulus