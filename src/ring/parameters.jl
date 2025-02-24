struct CyclicParam
    m::Int64
    N::Int64

    CyclicParam(m::Int64) = new(m, m)
end

struct CyclotomicParam
    m::Int64
    N::Int64

    CyclotomicParam(m::Int64) = begin
        @assert ispow2(m) || isodd(m) "$m is neither a power of two, nor an odd number."
        isprime(m) && @info "It is recommended to use the subring structure with d = 1 when the degree is a prime."
        new(m, totient(m))
    end
end

struct SubringParam
    m::Int64
    d::Int64
    N::Int64

    SubringParam(m::Int64, d::Int64) = begin
        @assert isprime(m) "$m is not a prime number."
        @assert (m - 1) % d == 0 "$d should divide $(m-1)."
        new(m, d, (m - 1) ÷ d)
    end
end

const RingParam = Union{CyclicParam,CyclotomicParam,SubringParam}

export CyclicParam, CyclotomicParam, SubringParam, RingParam

#==========================================================================#

function check_modulus(param::RingParam, Q::UInt64; onemod::Int64=1)
    if typeof(param) == CyclicParam
        _check_modulus_cyclic(param.m, Q; onemod=onemod)
    elseif typeof(param) == CyclotomicParam
        _check_modulus_cyclotomic(param.m, Q; onemod=onemod)
    elseif typeof(param) == SubringParam
        _check_modulus_subring(param.m, param.d, Q; onemod=onemod)
    end
end

function find_prime(param::RingParam, b::Real, n::Int64=1; onemod::Int64=1)
    if typeof(param) == CyclicParam
        _find_prime_cyclic(param.m, b, n; onemod=onemod)
    elseif typeof(param) == CyclotomicParam
        _find_prime_cyclotomic(param.m, b, n; onemod=onemod)
    elseif typeof(param) == SubringParam
        _find_prime_subring(param.m, param.d, b, n; onemod=onemod)
    end
end

function _check_modulus_cyclic(m::Int64, Q::UInt64; onemod::Int64)
    step = is2a3b5c7d(m) ? m : lcm(next2a3b5c7d(2m - 1), m)
    # Q % fact == 1
    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    _find_prime_cyclic(m::Int64, b::Real, n::Int64=1)

Return `n` primes of `b` bits, for cyclic NTT with degree `m`.
"""
function _find_prime_cyclic(m::Int64, b::Real, n::Int64=1; onemod::Int64=1)
    @assert b ≤ 62 "We do not support primes with bit length bigger than 62-bit"

    primes = Vector{UInt64}(undef, n)

    if is2a3b5c7d(m)
        step = lcm(m, onemod)
        initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
        currentn = nextprime(initialn, interval=step)
        initialn == currentn && (initialn -= step)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + step, interval=step)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=step)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - step, interval=step)
                cnt += 1
            end
        end
    else
        step = lcm(next2a3b5c7d(2m - 1), m, onemod)
        initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
        currentn = nextprime(initialn, interval=step)
        initialn == currentn && (initialn -= step)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + step, interval=step)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=step)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - step, interval=step)
                cnt += 1
            end
        end
    end

    primes
end

function _check_modulus_cyclotomic(m::Int64, Q::UInt64; onemod::Int64=1)
    if ispow2(m)
        step = lcm(m, onemod)
    else
        bluestein = next2a3b5c7d(2m - 1)

        N = totient(m)
        factors = factor(Vector, m)
        p = filter(isodd, factors)[1]
        moverp = m ÷ p
        deg = m - moverp
        if deg == N
            step = lcm(bluestein, m, onemod)
        else
            ñ = next2a3b5c7d(N)
            A = next2a3b5c7d(2(deg - N) + 1)
            step = lcm(bluestein, m, ñ, A, onemod)
        end
    end

    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    _find_prime_cyclotomic(m::Int64, b::Real, n::Int64=1)

Return `m` primes of `b`-bits, for cyclotomic NTT with degree `m`.
"""
function _find_prime_cyclotomic(m::Int64, b::Real, n::Int64=1; onemod::Int64=1)
    @assert b ≤ 62 "We do not support primes with bit length bigger than 62-bit"

    primes = Vector{UInt64}(undef, n)

    if ispow2(m)
        step = lcm(m, onemod)
        initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
        currentn = nextprime(initialn, interval=step)
        initialn == currentn && (initialn -= step)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + step, interval=step)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=step)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - step, interval=step)
                cnt += 1
            end
        end
    else
        bluestein = next2a3b5c7d(2m - 1)

        N = totient(m)
        factors = factor(Vector, m)
        p = filter(isodd, factors)[1]
        moverp = m ÷ p
        deg = m - moverp
        if deg == N
            step = lcm(bluestein, m, onemod)
        else
            ñ = next2a3b5c7d(N)
            A = next2a3b5c7d(2(deg - N) + 1)
            step = lcm(bluestein, m, ñ, A, onemod)
        end

        initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
        currentn = nextprime(initialn, interval=step)
        initialn == currentn && (initialn -= step)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + step, interval=step)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=step)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - step, interval=step)
                cnt += 1
            end
        end
    end

    primes
end

function _check_modulus_reductor(deg::Int64, N::Int64, Q::UInt64; onemod::Int64=1)
    ñ = next2a3b5c7d(N)
    A = next2a3b5c7d(2(deg - N) + 1)
    step = lcm(ñ, A, onemod)

    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    _find_prime_reductor(b::Real, deg::Int64, N::Int64, n::Int64=1)

Return `n` primes of `b`-bits, for reductor with degree N and at most degree deg.
"""
function _find_prime_reductor(deg::Int64, N::Int64, b::Real, n::Int64=1; onemod::Int64=1)
    @assert b ≤ 62 "We do not support primes with bit length bigger than 62-bit"

    primes = Vector{UInt64}(undef, n)

    ñ = next2a3b5c7d(N)
    A = next2a3b5c7d(2(deg - N) + 1)
    step = lcm(ñ, A, onemod)

    initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
    currentn = nextprime(initialn, interval=step)
    initialn == currentn && (initialn -= step)

    cnt = 1
    while log2(currentn) < 62 && cnt ≤ n
        primes[cnt] = currentn
        currentn = nextprime(currentn + step, interval=step)
        cnt += 1
    end

    if cnt ≤ n
        currentn = prevprime(initialn, interval=step)
        while cnt ≤ n
            primes[cnt] = currentn
            currentn = prevprime(currentn - step, interval=step)
            cnt += 1
        end
    end

    primes
end

function _check_modulus_subring(m::Int64, d::Int64, Q::UInt64; onemod::Int64=1)
    N = (m - 1) ÷ d
    step = is2a3b5c7d(N) ? lcm(N * m, onemod) : lcm(next2a3b5c7d(2N - 1), m, onemod)
    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    _find_prime_subring(m::Int64, p::Int64, b::Real, n::Int64=1)

Return `n` primes of `b`-bits for subring homomorphic encryption,
which allows the subring NTT with parameters `p` and `m`.
"""
function _find_prime_subring(m::Int64, d::Int64, b::Real, n::Int64=1; onemod::Int64=1)
    @assert b ≤ 62 "We do not support primes with bit length bigger than 62-bit"
    @assert isprime(m) "$m is not a prime number."

    N = (m - 1) ÷ d

    primes = Vector{UInt64}(undef, n)

    if is2a3b5c7d(N)
        step = lcm(N * m, onemod)
        initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
        currentn = nextprime(initialn, interval=step)
        initialn == currentn && (initialn -= step)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + step, interval=step)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=step)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - step, interval=step)
                cnt += 1
            end
        end
    else
        step = lcm(next2a3b5c7d(2N - 1), m, onemod)
        initialn = floor(UInt64, Double64(2.0)^b / step) * step + 1
        currentn = nextprime(initialn, interval=step)
        initialn == currentn && (initialn -= step)

        cnt = 1
        while log2(currentn) < 62 && cnt ≤ n
            primes[cnt] = currentn
            currentn = nextprime(currentn + step, interval=step)
            cnt += 1
        end

        if cnt ≤ n
            currentn = prevprime(initialn, interval=step)
            while cnt ≤ n
                primes[cnt] = currentn
                currentn = prevprime(currentn - step, interval=step)
                cnt += 1
            end
        end
    end

    primes
end

export check_modulus, find_prime