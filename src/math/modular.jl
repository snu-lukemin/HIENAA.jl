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
    halfQ::UInt64     # Used for balanced representation.

    function Modulus(Q::Integer)::Modulus
        Q == 0 && return new(zero(UInt64), zero(UInt64), zero(UInt64), zero(UInt64), zero(UInt64), 0, zero(UInt64))

        logQ = ceil(Int64, log2(Q))
        if logQ > 62
            throw(DomainError("Modulus should be smaller than 2⁶²."))
        end

        if isodd(Q)
            Q⁻¹ = UInt64(invmod(-Q, 0x00000000000000010000000000000000))
            R⁻¹ = UInt64(invmod(0x00000000000000010000000000000000, Q))
        else
            Q⁻¹ = zero(UInt64)
            R⁻¹ = zero(UInt64)
        end

        r = floor(UInt128, (big(1) << 128) / Q)
        r0, r1 = r % UInt64, (r >>> 64) % UInt64
        halfQ = (Q - 1) >>> 1
        new(Q, Q⁻¹, R⁻¹, r0, r1, logQ, halfQ)
    end
end

Base.:show(io::IO, Q::Modulus)::Nothing = begin
    print(io, Q.Q)
    return nothing
end
Base.:show(io::IO, mime::MIME"text/plain", Q::Modulus)::Nothing = begin
    println(io, Q.Q)
    return nothing
end

"""
Bred(x, Q) returns x % Q using barret reduction.
"""
Bred(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) >>> 64) % UInt64
    res = x - Q.Q * t
    res ≥ Q.Q ? res - Q.Q : res
end

Bred(x::Int64, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return Bred(UInt64(x), Q)
    else
        res = Bred(UInt64(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

Bred(x::UInt128, Q::Modulus)::UInt64 = begin
    x0 = x % UInt64
    x1 = (x >>> 64) % UInt64
    t = (widemul(Q.r1, x1) + (widemul(Q.r1, x0) + widemul(Q.r0, x1) + widemul(Q.r0, x0) >>> 64) >>> 64) % UInt64
    res = x0 - t * Q.Q
    res ≥ Q.Q ? res - Q.Q : res
end

Bred(x::Int128, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return Bred(UInt128(x), Q)
    else
        res = Bred(UInt128(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

lazy_Bred(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r1, x) >>> 64) % UInt64
    x - Q.Q * t
end

lazy_Bred(x::Int64, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return lazy_Bred(UInt64(x), Q)
    else
        res = lazy_Bred(UInt64(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

lazy_Bred(x::UInt128, Q::Modulus)::UInt64 = begin
    x0 = x % UInt64
    x1 = (x >>> 64) % UInt64
    t = (widemul(Q.r1, x1) + (widemul(Q.r1, x0) + widemul(Q.r0, x1) + widemul(Q.r0, x0) >>> 64) >>> 64) % UInt64
    x0 - t * Q.Q
end

lazy_Bred(x::Int128, Q::Modulus)::UInt64 = begin
    if x ≥ 0
        return lazy_Bred(UInt128(x), Q)
    else
        res = lazy_Bred(UInt128(-x), Q)
        return res > 0 ? Q.Q - res : res
    end
end

Bred!(a::AbstractVector{UInt64}, Q::Modulus)::Nothing = Bred_to!(a, a, Q)

Bred_to!(res::AbstractVector{UInt64}, a::AbstractVector{<:Union{Int64,UInt64,Int128,UInt128}}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = Bred(a[i+1], Q)
        res[i+2] = Bred(a[i+2], Q)
        res[i+3] = Bred(a[i+3], Q)
        res[i+4] = Bred(a[i+4], Q)
        res[i+5] = Bred(a[i+5], Q)
        res[i+6] = Bred(a[i+6], Q)
        res[i+7] = Bred(a[i+7], Q)
        res[i+8] = Bred(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bred(a[i], Q)
    end

    return nothing
end

# Taken from the invmod.
Base.:invmod(x::Integer, Q::Modulus)::UInt64 = begin
    n, m = UInt64(x), Q.Q
    g, x, _ = gcdx(n, m)
    g ≠ 1 && throw(DomainError((n, m), LazyString("Greatest common divisor is ", g, ".")))
    x > typemax(UInt64) >>> 1 && (x += m)
    Bred(x, Q)
end

Mform(x::UInt64, Q::Modulus)::UInt64 = begin
    t = (widemul(Q.r0, x) >>> 64) % UInt64 + Q.r1 * x
    res = -t * Q.Q
    res ≥ Q.Q ? res - Q.Q : res
end

Mform(x::Int64, Q::Modulus)::UInt64 = Mform(Bred(x, Q), Q)

iMform(x::UInt64, Q::Modulus)::UInt64 = Bred(widemul(x, Q.R⁻¹), Q)

Mform!(x::AbstractVector{UInt64}, Q::Modulus)::Nothing = Mform_to!(x, x, Q)

Mform_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = Mform(a[i+1], Q)
        res[i+2] = Mform(a[i+2], Q)
        res[i+3] = Mform(a[i+3], Q)
        res[i+4] = Mform(a[i+4], Q)
        res[i+5] = Mform(a[i+5], Q)
        res[i+6] = Mform(a[i+6], Q)
        res[i+7] = Mform(a[i+7], Q)
        res[i+8] = Mform(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Mform(a[i], Q)
    end

    return nothing
end

iMform!(x::AbstractVector{UInt64}, Q::Modulus)::Nothing = iMform_to!(x, x, Q)

iMform_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = iMform(a[i+1], Q)
        res[i+2] = iMform(a[i+2], Q)
        res[i+3] = iMform(a[i+3], Q)
        res[i+4] = iMform(a[i+4], Q)
        res[i+5] = iMform(a[i+5], Q)
        res[i+6] = iMform(a[i+6], Q)
        res[i+7] = iMform(a[i+7], Q)
        res[i+8] = iMform(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = iMform(a[i], Q)
    end

    return nothing
end

"""
Convert UInt64 to a balanced representation.
"""
Base.:signed(x::UInt64, Q::Modulus)::Int64 = signed(x > Q.halfQ ? x - Q.Q : x)

add(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    res = a + b
    res ≥ Q.Q ? res - Q.Q : res
end

add_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = add(a[i+1], b[i+1], Q)
        res[i+2] = add(a[i+2], b[i+2], Q)
        res[i+3] = add(a[i+3], b[i+3], Q)
        res[i+4] = add(a[i+4], b[i+4], Q)
        res[i+5] = add(a[i+5], b[i+5], Q)
        res[i+6] = add(a[i+6], b[i+6], Q)
        res[i+7] = add(a[i+7], b[i+7], Q)
        res[i+8] = add(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = add(a[i], b[i], Q)
    end

    return nothing
end

add_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::UInt64, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = add(a[i+1], b, Q)
        res[i+2] = add(a[i+2], b, Q)
        res[i+3] = add(a[i+3], b, Q)
        res[i+4] = add(a[i+4], b, Q)
        res[i+5] = add(a[i+5], b, Q)
        res[i+6] = add(a[i+6], b, Q)
        res[i+7] = add(a[i+7], b, Q)
        res[i+8] = add(a[i+8], b, Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = add(a[i], b, Q)
    end

    return nothing
end

neg(a::UInt64, Q::Modulus)::UInt64 = a == 0 ? a : Q.Q - a

neg_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = neg(a[i+1], Q)
        res[i+2] = neg(a[i+2], Q)
        res[i+3] = neg(a[i+3], Q)
        res[i+4] = neg(a[i+4], Q)
        res[i+5] = neg(a[i+5], Q)
        res[i+6] = neg(a[i+6], Q)
        res[i+7] = neg(a[i+7], Q)
        res[i+8] = neg(a[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = neg(a[i], Q)
    end

    return nothing
end

sub(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    res = a - b
    res ≥ Q.Q ? res + Q.Q : res
end

sub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = sub(a[i+1], b[i+1], Q)
        res[i+2] = sub(a[i+2], b[i+2], Q)
        res[i+3] = sub(a[i+3], b[i+3], Q)
        res[i+4] = sub(a[i+4], b[i+4], Q)
        res[i+5] = sub(a[i+5], b[i+5], Q)
        res[i+6] = sub(a[i+6], b[i+6], Q)
        res[i+7] = sub(a[i+7], b[i+7], Q)
        res[i+8] = sub(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = sub(a[i], b[i], Q)
    end

    return nothing
end

sub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::UInt64, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = sub(a[i+1], b, Q)
        res[i+2] = sub(a[i+2], b, Q)
        res[i+3] = sub(a[i+3], b, Q)
        res[i+4] = sub(a[i+4], b, Q)
        res[i+5] = sub(a[i+5], b, Q)
        res[i+6] = sub(a[i+6], b, Q)
        res[i+7] = sub(a[i+7], b, Q)
        res[i+8] = sub(a[i+8], b, Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = sub(a[i], b, Q)
    end

    return nothing
end

Mmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    q, q⁻¹ = Q.Q, Q.Q⁻¹

    ab = widemul(a, b)
    w = ((widemul(q, (ab % UInt64) * q⁻¹) + ab) >>> 64) % UInt64
    w ≥ q ? w - q : w
end

lazy_Mmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = begin
    q, q⁻¹ = Q.Q, Q.Q⁻¹

    ab = widemul(a, b)
    ((widemul(q, (ab % UInt64) * q⁻¹) + ab) >>> 64) % UInt64
end

Bmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = Bred(widemul(a, b), Q)

lazy_Bmul(a::UInt64, b::UInt64, Q::Modulus)::UInt64 = lazy_Bred(widemul(a, b), Q)

Mmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = Mmul(a[i+1], b[i+1], Q)
        res[i+2] = Mmul(a[i+2], b[i+2], Q)
        res[i+3] = Mmul(a[i+3], b[i+3], Q)
        res[i+4] = Mmul(a[i+4], b[i+4], Q)
        res[i+5] = Mmul(a[i+5], b[i+5], Q)
        res[i+6] = Mmul(a[i+6], b[i+6], Q)
        res[i+7] = Mmul(a[i+7], b[i+7], Q)
        res[i+8] = Mmul(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Mmul(a[i], b[i], Q)
    end 

    return nothing
end

Mmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = Mmul(a, b[i+1], Q)
        res[i+2] = Mmul(a, b[i+2], Q)
        res[i+3] = Mmul(a, b[i+3], Q)
        res[i+4] = Mmul(a, b[i+4], Q)
        res[i+5] = Mmul(a, b[i+5], Q)
        res[i+6] = Mmul(a, b[i+6], Q)
        res[i+7] = Mmul(a, b[i+7], Q)
        res[i+8] = Mmul(a, b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Mmul(a, b[i], Q)
    end

    return nothing
end

Bmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = Bmul(a[i+1], b[i+1], Q)
        res[i+2] = Bmul(a[i+2], b[i+2], Q)
        res[i+3] = Bmul(a[i+3], b[i+3], Q)
        res[i+4] = Bmul(a[i+4], b[i+4], Q)
        res[i+5] = Bmul(a[i+5], b[i+5], Q)
        res[i+6] = Bmul(a[i+6], b[i+6], Q)
        res[i+7] = Bmul(a[i+7], b[i+7], Q)
        res[i+8] = Bmul(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bmul(a[i], b[i], Q)
    end

    return nothing
end

Bmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = Bmul(a, b[i+1], Q)
        res[i+2] = Bmul(a, b[i+2], Q)
        res[i+3] = Bmul(a, b[i+3], Q)
        res[i+4] = Bmul(a, b[i+4], Q)
        res[i+5] = Bmul(a, b[i+5], Q)
        res[i+6] = Bmul(a, b[i+6], Q)
        res[i+7] = Bmul(a, b[i+7], Q)
        res[i+8] = Bmul(a, b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bmul(a, b[i], Q)
    end

    return nothing
end

lazy_Bmul_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = lazy_Bmul(a[i+1], b[i+1], Q)
        res[i+2] = lazy_Bmul(a[i+2], b[i+2], Q)
        res[i+3] = lazy_Bmul(a[i+3], b[i+3], Q)
        res[i+4] = lazy_Bmul(a[i+4], b[i+4], Q)
        res[i+5] = lazy_Bmul(a[i+5], b[i+5], Q)
        res[i+6] = lazy_Bmul(a[i+6], b[i+6], Q)
        res[i+7] = lazy_Bmul(a[i+7], b[i+7], Q)
        res[i+8] = lazy_Bmul(a[i+8], b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = lazy_Bmul(a[i], b[i], Q)
    end

    return nothing
end

lazy_Bmul_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = lazy_Bmul(a, b[i+1], Q)
        res[i+2] = lazy_Bmul(a, b[i+2], Q)
        res[i+3] = lazy_Bmul(a, b[i+3], Q)
        res[i+4] = lazy_Bmul(a, b[i+4], Q)
        res[i+5] = lazy_Bmul(a, b[i+5], Q)
        res[i+6] = lazy_Bmul(a, b[i+6], Q)
        res[i+7] = lazy_Bmul(a, b[i+7], Q)
        res[i+8] = lazy_Bmul(a, b[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = lazy_Bmul(a, b[i], Q)
    end

    return nothing
end

Bmuladd_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = Bred(widemul(a[i+1], b[i+1]) + res[i+1], Q)
        res[i+2] = Bred(widemul(a[i+2], b[i+2]) + res[i+2], Q)
        res[i+3] = Bred(widemul(a[i+3], b[i+3]) + res[i+3], Q)
        res[i+4] = Bred(widemul(a[i+4], b[i+4]) + res[i+4], Q)
        res[i+5] = Bred(widemul(a[i+5], b[i+5]) + res[i+5], Q)
        res[i+6] = Bred(widemul(a[i+6], b[i+6]) + res[i+6], Q)
        res[i+7] = Bred(widemul(a[i+7], b[i+7]) + res[i+7], Q)
        res[i+8] = Bred(widemul(a[i+8], b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bred(widemul(a[i], b[i]) + res[i], Q)
    end

    return nothing
end

Bmuladd_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(b)
    @inbounds for i = 0:8:N-8
        res[i+1] = Bred(widemul(a, b[i+1]) + res[i+1], Q)
        res[i+2] = Bred(widemul(a, b[i+2]) + res[i+2], Q)
        res[i+3] = Bred(widemul(a, b[i+3]) + res[i+3], Q)
        res[i+4] = Bred(widemul(a, b[i+4]) + res[i+4], Q)
        res[i+5] = Bred(widemul(a, b[i+5]) + res[i+5], Q)
        res[i+6] = Bred(widemul(a, b[i+6]) + res[i+6], Q)
        res[i+7] = Bred(widemul(a, b[i+7]) + res[i+7], Q)
        res[i+8] = Bred(widemul(a, b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bred(widemul(a, b[i]) + res[i], Q)
    end

    return nothing
end

Bmulsub_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(a) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(a)
    @inbounds for i = 0:8:N-8
        res[i+1] = Bred(widemul(Q.Q - a[i+1], b[i+1]) + res[i+1], Q)
        res[i+2] = Bred(widemul(Q.Q - a[i+2], b[i+2]) + res[i+2], Q)
        res[i+3] = Bred(widemul(Q.Q - a[i+3], b[i+3]) + res[i+3], Q)
        res[i+4] = Bred(widemul(Q.Q - a[i+4], b[i+4]) + res[i+4], Q)
        res[i+5] = Bred(widemul(Q.Q - a[i+5], b[i+5]) + res[i+5], Q)
        res[i+6] = Bred(widemul(Q.Q - a[i+6], b[i+6]) + res[i+6], Q)
        res[i+7] = Bred(widemul(Q.Q - a[i+7], b[i+7]) + res[i+7], Q)
        res[i+8] = Bred(widemul(Q.Q - a[i+8], b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bred(widemul(Q.Q - a[i], b[i]) + res[i], Q)
    end

    return nothing
end

Bmulsub_to!(res::AbstractVector{UInt64}, a::UInt64, b::AbstractVector{UInt64}, Q::Modulus)::Nothing = begin
    if length(res) ≠ length(b)
        throw(DimensionMismatch("Lengths of input and output vectors should be the same."))
    end
    N = length(b)
    neg_a = Q.Q - a
    @inbounds for i = 0:8:N-8
        res[i+1] = Bred(widemul(neg_a, b[i+1]) + res[i+1], Q)
        res[i+2] = Bred(widemul(neg_a, b[i+2]) + res[i+2], Q)
        res[i+3] = Bred(widemul(neg_a, b[i+3]) + res[i+3], Q)
        res[i+4] = Bred(widemul(neg_a, b[i+4]) + res[i+4], Q)
        res[i+5] = Bred(widemul(neg_a, b[i+5]) + res[i+5], Q)
        res[i+6] = Bred(widemul(neg_a, b[i+6]) + res[i+6], Q)
        res[i+7] = Bred(widemul(neg_a, b[i+7]) + res[i+7], Q)
        res[i+8] = Bred(widemul(neg_a, b[i+8]) + res[i+8], Q)
    end
    @inbounds for i = 8(N>>3)+1:N
        res[i] = Bred(widemul(neg_a, b[i]) + res[i], Q)
    end

    return nothing
end

const Moduli::DataType = AbstractVector{Modulus}

Base.:prod(Q::Moduli)::BigInt = begin
    Qbig = BigInt(1)
    @inbounds for i = eachindex(Q)
        Qbig *= Q[i].Q
    end
    Qbig
end

Base.:log2(Q::Moduli)::Float64 = begin
    logQ = 0.0
    @inbounds for i = eachindex(Q)
        logQ += log2(Q[i].Q)
    end
    logQ
end