abstract type RingParam end

struct CyclicParam <: RingParam
    m::Int64
    N::Int64

    CyclicParam(m::Int64)::CyclicParam = new(m, m)
end

struct CyclotomicParam <: RingParam
    m::Int64
    N::Int64

    CyclotomicParam(m::Int64)::CyclotomicParam = begin
        if !ispow2(m) && !isodd(m)
            throw(DomainError("Cyclotomic degree should be either a power of two, or an odd number."))
        end
        isprime(m) && @info "It is recommended to use the subring structure with d = 1 when the degree is a prime."
        new(m, totient(m))
    end
end

struct SubringParam <: RingParam
    m::Int64
    d::Int64
    N::Int64

    SubringParam(m::Int64, d::Int64)::SubringParam = begin
        if !isprime(m)
            throw(DomainError("The degree should be a prime number."))
        end
        if (m - 1) % d ≠ 0
            throw(DomainError("d should divide $(m-1)."))
        end
        new(m, d, (m - 1) ÷ d)
    end
end

#==========================================================================#

function check_modulus(param::RingParam, Q::UInt64; onemod::Int64=1)::Bool
    if isa(param, CyclicParam)
        check_modulus_cyclic(param.m, Q; onemod=onemod)
    elseif isa(param, CyclotomicParam)
        check_modulus_cyclotomic(param.m, Q; onemod=onemod)
    elseif isa(param, SubringParam)
        check_modulussubring(param.m, param.d, Q; onemod=onemod)
    end
end

function find_prime(param::RingParam, b::Real, n::Int64=1; onemod::Int64=1)::Vector{UInt64}
    if isa(param, CyclicParam)
        find_prime_cyclic(param.m, b, n; onemod=onemod)
    elseif isa(param, CyclotomicParam)
        find_prime_cyclotomic(param.m, b, n; onemod=onemod)
    elseif isa(param, SubringParam)
        find_primesubring(param.m, param.d, b, n; onemod=onemod)
    end
end

function check_modulus_cyclic(m::Int64, Q::UInt64; onemod::Int64)::Bool
    step = is2a3b5c7d(m) ? m : lcm(next2a3b(2m - 1), m)
    # Q % fact == 1
    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    find_prime_cyclic(m::Int64, b::Real, n::Int64=1)

Return `n` primes of `b` bits, for cyclic NTT with degree `m`.
"""
function find_prime_cyclic(m::Int64, b::Real, n::Int64=1; onemod::Int64=1)::Vector{UInt64}
    if b > 62
        throw(DomainError("We do not support primes with bit length bigger than 62-bit"))
    end

    primes = Vector{UInt64}(undef, n)

    if is2a3b5c7d(m)
        step = lcm(m, onemod)
        initialn = floor(UInt64, 2.0^b / step) * step + 1
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
        step = lcm(next2a3b(2m - 1), m, onemod)
        initialn = floor(UInt64, 2.0^b / step) * step + 1
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

function check_modulus_cyclotomic(m::Int64, Q::UInt64; onemod::Int64=1)::Bool
    if ispow2(m)
        step = lcm(m, onemod)
    else
        bluestein = next2a3b(2m - 1)

        N = totient(m)
        factors = factor(Vector, m)
        p = filter(isodd, factors)[1]
        moverp = m ÷ p
        deg = m - moverp
        if deg == N
            step = lcm(bluestein, m, onemod)
        else
            ñ = next2a3b(N)
            A = next2a3b(2(deg - N) + 1)
            step = lcm(bluestein, m, ñ, A, onemod)
        end
    end

    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    find_prime_cyclotomic(m::Int64, b::Real, n::Int64=1)

Return `m` primes of `b`-bits, for cyclotomic NTT with degree `m`.
"""
function find_prime_cyclotomic(m::Int64, b::Real, n::Int64=1; onemod::Int64=1)::Vector{UInt64}
    if b > 62
        throw(DomainError("We do not support primes with bit length bigger than 62-bit"))
    end

    primes = Vector{UInt64}(undef, n)

    if ispow2(m)
        step = lcm(m, onemod)
        initialn = floor(UInt64, 2.0^b / step) * step + 1
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
        bluestein = next2a3b(2m - 1)

        N = totient(m)
        factors = factor(Vector, m)
        p = filter(isodd, factors)[1]
        moverp = m ÷ p
        deg = m - moverp
        if deg == N
            step = lcm(bluestein, m, onemod)
        else
            ñ = next2a3b(N)
            A = next2a3b(2(deg - N) + 1)
            step = lcm(bluestein, m, ñ, A, onemod)
        end

        initialn = floor(UInt64, 2.0^b / step) * step + 1
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

function check_modulus_reductor(deg::Int64, N::Int64, Q::UInt64; onemod::Int64=1)::Bool
    ñ = next2a3b(N)
    A = next2a3b(2(deg - N) + 1)
    step = lcm(ñ, A, onemod)

    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    find_prime_reductor(b::Real, deg::Int64, N::Int64, n::Int64=1)

Return `n` primes of `b`-bits, for reductor with degree N and at most degree deg.
"""
function find_prime_reductor(deg::Int64, N::Int64, b::Real, n::Int64=1; onemod::Int64=1)::Vector{UInt64}
    if b > 62
        throw(DomainError("We do not support primes with bit length bigger than 62-bit"))
    end

    primes = Vector{UInt64}(undef, n)

    ñ = next2a3b(N)
    A = next2a3b(2(deg - N) + 1)
    step = lcm(ñ, A, onemod)

    initialn = floor(UInt64, 2.0^b / step) * step + 1
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

function check_modulussubring(m::Int64, d::Int64, Q::UInt64; onemod::Int64=1)::Bool
    N = (m - 1) ÷ d
    step = is2a3b5c7d(N) ? lcm(N * m, onemod) : lcm(next2a3b(2N - 1), m, onemod)
    Qi = collect(keys(factor(Dict, Q)))
    all(Qi .% step .== 1)
end

"""
    find_primesubring(m::Int64, p::Int64, b::Real, n::Int64=1)

Return `n` primes of `b`-bits for subring homomorphic encryption,
which allows the subring NTT with parameters `p` and `m`.
"""
function find_primesubring(m::Int64, d::Int64, b::Real, n::Int64=1; onemod::Int64=1)::Vector{UInt64}
    if b > 62
        throw(DomainError("We do not support primes with bit length bigger than 62-bit"))
    end
    if !isprime(m)
        throw(DomainError("The degree should be a prime number."))
    end

    N = (m - 1) ÷ d

    primes = Vector{UInt64}(undef, n)

    if is2a3b5c7d(N)
        step = lcm(N * m, onemod)
        initialn = floor(UInt64, 2.0^b / step) * step + 1
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
        step = lcm(next2a3b(2N - 1), m, onemod)
        initialn = floor(UInt64, 2.0^b / step) * step + 1
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