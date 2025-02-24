"""
    ord(p::Integer, m::Integer)

Returns the multiplicative order of p modulo m.
"""
function ord(p::Integer, m::Integer)::Int64
    res = 1
    acc = p % m
    while acc ≠ 1
        acc = (acc * p) % m
        res += 1
    end

    res
end

zeropadto(a::Vector{Int64}, n::Int64) = begin
    @assert n ≥ length(a)
    vcat(a, zeros(Int64, n - length(a)))
end

zeropadto(a::Vector{UInt64}, n::Int64) = begin
    @assert n ≥ length(a)
    vcat(a, zeros(UInt64, n - length(a)))
end

@inline ispow3(n::Int64)::Bool = n == 3^round(Int64, log(n) / log(3))

@inline ispow5(n::Int64)::Bool = n == 5^round(Int64, log(n) / log(5))

@inline ispow7(n::Int64)::Bool = n == 7^round(Int64, log(n) / log(7))

@noinline is2a3b5c7d(n::Int64)::Bool = keys(factor(Dict, n)) ⊆ [2, 3, 5, 7]

@noinline next2a3b5c7d(n::Int64)::Int64 = begin
    res = n + 1
    while !is2a3b5c7d(res)
        res += 1
    end
    res
end

@inline factor2357(n::Int64)::NTuple{4,Int64} = begin
    n2, n3, n5, n7 = 1, 1, 1, 1

    while n % 2 == 0
        n >>= 1
        n2 <<= 1
    end

    while n % 3 == 0
        n ÷= 3
        n3 *= 3
    end

    while n % 5 == 0
        n ÷= 5
        n5 *= 5
    end

    while n % 7 == 0
        n ÷= 7
        n7 *= 7
    end

    n2, n3, n5, n7
end

# This parameter guarantees at most -52 b of fixed point error (with 64 levels). 
# The decryption failure will happen in a very unlikely situation...
const fixed_prec = 124
const float_prec = 62
const fixed_mask = UInt64(1) << float_prec - 1
const round_mask = UInt64(1) << (fixed_prec - float_prec) - 1

mult_and_round(a::UInt64, b::UInt128) = a * (b >> float_prec) + ((a * (b & fixed_mask)) >> float_prec)
round_to_uint64(a::UInt128) = round_to_uint128(a) % UInt64
round_to_uint128(a::UInt128) = (a >> (fixed_prec - float_prec) + (a & round_mask) >> (fixed_prec - float_prec - 1))
round_to_uint128(a::Int128) = unsigned((a >> (fixed_prec - float_prec) + (a & round_mask) >> (fixed_prec - float_prec - 1)))

"""
scramble!(a, r) computes an in-place r-radix reversal algorithm. 
The length of the input vector should be a power-of-r.
"""
function scramble!(a::AbstractVector{<:Number}, r::Int64)
    if r == 2
        j = 0
        N = length(a)
        @inbounds for i = 1:N-1
            bit = N >> 1
            while j ≥ bit
                j -= bit
                bit >>= 1
            end
            j += bit
            if i < j
                a[i+1], a[j+1] = a[j+1], a[i+1]
            end
        end
    else
        logN = round(Int64, log(length(a)) / log(r))
        N = r^logN

        @inbounds for i = 1:N-1
            idx, j = i, 0
            @simd for _ = 1:logN
                j = j * r + (idx % r)
                idx ÷= r
            end
            if j > i
                a[i+1], a[j+1] = a[j+1], a[i+1]
            end
        end
    end
end

Base.:sum(a::AbstractVector{UInt64}, Q::Modulus)::UInt64 = begin
    res = zero(UInt64)
    @inbounds for i = eachindex(a)
        res = _add(res, a[i], Q)
    end
    res
end

Base.:powermod(a::UInt64, p::Integer, Q::Modulus)::UInt64 = begin
    @assert p ≥ 0
    p == 0 && return UInt64(1)

    t = prevpow(2, p)
    r = UInt64(1)
    while true
        if p >= t
            r = _Bmul(r, a, Q)
            p -= t
        end
        t >>>= 1
        t <= 0 && break
        r = _Bmul(r, r, Q)
    end
    return r
end

"""
primitive_root_finder(Q) finds the primitive root of Q.
"""
function primitive_root_finder(Q::Modulus)
    facts = factor(Dict, Q.Q)
    Qi = keys(facts) .^ values(facts)

    if length(Qi) == 1
        N = totient(Q.Q)
        test = N .÷ collect(keys(factor(Dict, N)))

        g = UInt64(2)
        while true
            if powermod(g, N, Q) == 1 && all(powermod.(g, test, Ref(Q)) .≠ 1)
                break
            end
            g += UInt64(1)
        end

        g
    else
        gi = primitive_root_finder.(Qi)

        res = UInt64(0)
        for i = eachindex(Qi)
            tildei = Q.Q ÷ Qi[i]
            invtilde = invmod(tildei, Qi[i])

            res = _add(res, _Bmul(_Bmul(gi[i], invtilde, Q), tildei, Q), Q)
        end

        res
    end
end

primitive_root_finder(Q::Integer)::UInt64 = primitive_root_finder(Modulus(Q))

function ith_root_from_primitive_root(i::Int64, ξ::UInt64, Q::Modulus)
    facts = factor(Dict, Q.Q)
    Qiq = collect(keys(facts))
    Qir = collect(values(facts))
    Qi = Qiq .^ Qir

    ζ = UInt64(0)
    for idx = eachindex(Qi)
        ξi = ξ % Qi[idx]
        ζi = powermod(ξi, (Qiq[idx]^Qir[idx] - Qiq[idx]^(Qir[idx] - 1)) ÷ i, Qi[idx])
        tildei = Q.Q ÷ Qi[idx]
        invtilde = invmod(tildei, Qi[idx])

        ζ = _add(ζ, _Bmul(_Bmul(ζi, invtilde, Q), tildei, Q), Q)
    end

    ζ
end

function find_generators_mod_m(m::Int64)
    fact_dict = factor(Dict, m)
    facs = collect(fact_dict)
    sort!(facs, by=x -> x[1])

    dims = Memory{Int64}(undef, length(facs))
    gens = Memory{Int64}(zeros(Int64, length(facs)))

    @inbounds for i = eachindex(facs)
        p, e = facs[i].first, facs[i].second

        for j = eachindex(facs)
            i == j && continue
            pj, ej = facs[j].first, facs[j].second
            gens[i] += invmod(m ÷ pj^ej, pj^ej) * m ÷ pj^ej
        end

        dims[i] = p^e - p^(e - 1)

        for j = 2:p^e-1
            j % p == 0 && continue
            if ord(j, p^e) == dims[i]
                gens[i] += j * invmod(m ÷ p^e, p^e) * m ÷ p^e
                break
            end
        end

        gens[i] %= m
    end

    dims, gens
end

"""
division(a, b) performs the long division of polynomials.
This function affects the vector a, so be cautious when using.
"""
function division(a::Vector{Int64}, b::Vector{Int64})::Vector{Int64}
    @assert length(a) ≥ length(b)

    Q = zeros(Int64, length(a) - length(b) + 1)

    @inbounds for i = 0:length(a)-length(b)
        if a[end-i] ≠ 0
            Q[end-i] = a[end-i] ÷ b[end]

            @simd for j = 0:length(b)-1
                a[end-i-j] -= b[end-j] * Q[end-i]
            end
        end
    end

    Q
end

"""
division(a, b) performs the long division of polynomials.
This function affects the vector a, so be cautious when using.
"""
function division_mod_Q(a::Vector{UInt64}, b::Vector{UInt64}, Q::Modulus)
    @assert length(a) ≥ length(b) "The size of the input polynomials do not match."
    @assert gcd(b[end], Q.Q) == 1 "The leading coefficient of b should be coprime with Q."

    q = zeros(UInt64, length(a) - length(b) + 1)

    @inbounds for i = 0:length(a)-length(b)
        if a[end-i] ≠ 0
            q[end-i] = _Bmul(a[end-i], invmod(b[end], Q.Q), Q)

            @simd for j = 0:length(b)-1
                a[end-i-j] = _neg(_Bred(widemul(b[end-j], q[end-i]) + Q.Q - a[end-i-j], Q), Q)
            end
        end
    end

    q
end

"""
cyclotomic_finder(m) returns the m-th cyclotomic polynomial. 
"""
function cyclotomic_finder(m::Int64)
    # Find the factors of m. Sorting is necessary to enhance the efficiency.
    factors = factor(Dict, m)
    primes = collect(keys(factors))
    sort!(primes)

    if primes[1] == 2
        is2divm = true
        popfirst!(primes)
    else
        is2divm = false
    end

    # Make arrays for less allocations.
    ϕm = zeros(Int64, m + 1)
    ϕm[1] = -1
    ϕm[2] = 1
    tmp = Vector{Int64}(undef, m + 1)
    tmp2 = Vector{Int64}(undef, m + 1)

    # Perform naive polynomial long divisions.
    η = 1
    ηprime = 1
    @inbounds for p = primes
        @. tmp = ϕm
        @. tmp2 = 0
        @. ϕm = 0
        @simd for i = 0:η
            tmp2[i*p+1] = tmp[i+1]
        end

        ηtimesp = η * p
        ηprime = ηtimesp - η

        for i = 0:(p-1)*η
            if tmp2[ηtimesp+1-i] ≠ 0
                ϕm[ηprime+1-i] = tmp2[ηtimesp+1-i] ÷ tmp[η+1]

                @simd for j = 0:η
                    tmp2[ηtimesp+1-i-j] -= tmp[η+1-j] * ϕm[ηprime+1-i]
                end
            end
        end

        η = ηprime
    end

    s = m ÷ reduce(*, primes)

    if is2divm
        @inbounds @simd for i = 1:2:η
            ϕm[i+1] = -ϕm[i+1]
        end

        s >>= 1
    end

    if s > 1
        @. tmp = ϕm
        @. ϕm = 0
        @inbounds @simd for i = 0:η
            ϕm[s*i+1] = tmp[i+1]
        end
    end

    resize!(ϕm, totient(m) + 1)
    ϕm
end

function gen_power_modmul(gen::Vector{Int64}, dims::Vector{Int64}, idx::NTuple{N,Int64}, m::Int64) where {N}
    if N == 1
        idx1 = powermod(gen[1], idx[1] % dims[1], m)
        res = idx1
    elseif N == 2
        idx1 = powermod(gen[1], idx[1] % dims[1], m)
        idx2 = powermod(gen[2], idx[2] % dims[2], m)
        res = (idx1 * idx2) % m
    elseif N == 3
        idx1 = powermod(gen[1], idx[1] % dims[1], m)
        idx2 = powermod(gen[2], idx[2] % dims[2], m)
        idx3 = powermod(gen[3], idx[3] % dims[3], m)
        res = (((idx1 * idx2) % m) * idx3) % m
    else
        idx1 = powermod(gen[1], idx[1] % dims[1], m)
        idx2 = powermod(gen[2], idx[2] % dims[2], m)
        idx3 = powermod(gen[3], idx[3] % dims[3], m)
        idx4 = powermod(gen[4], idx[4] % dims[4], m)
        res = (((((idx1 * idx2) % m) * idx3) % m) * idx4) % m
    end

    res
end