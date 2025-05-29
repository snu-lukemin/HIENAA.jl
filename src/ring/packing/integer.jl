"""
    IntPacker(t::Integer, param::RingParam)

IntPacker is an abstract type which supports SIMD packing modulo t.
"""
abstract type IntPacker end

struct IntPackerPow2 <: IntPacker
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    pr::Modulus
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    buff::Vector{UInt64}
    ntter::CyclotomicNTTransformerPow2

    function IntPackerPow2(m::Integer, p::Integer, r::Integer)::IntPackerPow2
        if !isprime(p) || iseven(p)
            throw(DomainError("Only odd prime power plaintext moduli are supported."))
        end

        d, N = ord(p, m), m >> 1
        k = N ÷ d
        pr = Modulus(p^r)

        if p % 8 == 1
            cube = [k >> 1, 2]
            cubegen = [5, 2N - 1]
            buff = Vector{UInt64}(undef, k)
            ntter = CyclotomicNTTransformerPow2(k << 1, pr)

            new(m, N, d, k, pr, cube, cubegen, buff, ntter)
        else
            error("Only the primes of the form 8k+1 are supported now.")
        end
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, msg::AbstractVector{UInt64}, packer::IntPackerPow2)::Nothing
    if length(res) ≠ packer.N
        throw(DimensionMismatch("The length of the output message should be equal to the polynomial degree."))
    end
    if length(msg) ≠ packer.k
        throw(DimensionMismatch("The length of the input message should be equal to the number of slots."))
    end

    k, d, ntter = packer.k, packer.d, packer.ntter

    @. res = 0
    tmp, mask = 1, 2k - 1
    @inbounds for i = 0:k>>1-1
        res[(tmp>>1)*d+1] = msg[k>>1-i]
        res[(k-tmp>>1-1)*d+1] = msg[k-i]
        tmp = 5tmp & mask
    end
    scramble!(res[1:d:end], 2)
    intt!(res[1:d:end], ntter)

    return nothing
end

@views function unpack_to!(res::AbstractVector{UInt64}, pt::AbstractVector{UInt64}, packer::IntPackerPow2)::Nothing
    if length(res) ≠ packer.k
        throw(DimensionMismatch("The length of the output message should be equal to the number of slots."))
    end
    if length(pt) ≠ packer.N
        throw(DimensionMismatch("The length of the input message should be equal to the polynomial degree."))
    end

    k, d, buff, ntter = packer.k, packer.d, packer.buff, packer.ntter

    @. buff = pt[1:d:end]
    ntt!(buff, ntter)
    scramble!(buff, 2)
    tmp, mask = 1, 2k - 1
    @inbounds for i = 0:k>>1-1
        res[k>>1-i] = buff[(tmp>>1)+1]
        res[k-i] = buff[(k-tmp>>1-1)+1]
        tmp = 5tmp & mask
    end

    return nothing
end

struct IntPackerNTT <: IntPacker
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    pr::Modulus
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    ntter::CyclotomicNTTransformerArb

    function IntPackerNTT(m::Integer, p::Integer, r::Integer)::IntPackerNTT
        if !isprime(p)
            throw(DomainError("Only prime power plaintext moduli are supported."))
        end
        N = totient(m)
        pr = Modulus(p^r)
        ntter = CyclotomicNTTransformerArb(m, pr)

        new(m, N, 1, N, pr, ntter.dims, ntter.gens, ntter)
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, msg::AbstractVector{UInt64}, packer::IntPackerNTT)::Nothing
    if !(length(res) == length(msg) == packer.N)
        throw(DimensionMismatch("The length of the input and output messages should be equal to the polynomial degree."))
    end

    @. res = msg
    intt!(res, packer.ntter)

    return nothing
end

@views function unpack_to!(res::AbstractVector{UInt64}, pt::AbstractVector{UInt64}, packer::IntPackerNTT)::Nothing
    if !(length(res) == length(pt) == packer.N)
        throw(DimensionMismatch("The length of the input and output messages should be equal to the polynomial degree."))
    end

    @. res = pt
    ntt!(res, packer.ntter)

    return nothing
end

# Many thanks to Simon Pohmann for the helpful discussions.
# We determine slot structure using an algorithm from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8691744
function factor_cyclotomic(m::Int64, p::Int64, r::Int64)::Tuple{Vector{Int64},Vector{Int64},Vector{Vector{Int64}},Vector{UInt64}}
    N = totient(m)
    d = ord(p, m)
    k = N ÷ d

    Kp, _ = finite_field(p, d, " ")
    ZZy, _ = ZZ["y"]
    ZZpr, _ = residue_ring(ZZ, p^r)
    ZZprz, z = polynomial_ring(ZZpr, "z")

    poly = ZZprz(lift(ZZy, defining_polynomial(Kp)))
    Kpr, _ = residue_ring(ZZprz, poly)
    Kprx, x = Kpr["x"]

    # Order of the multiplicative group of the multiplicative ring.
    order = big(p)^((r - 1) * d) * (big(p)^d - 1)
    test = m .÷ collect(keys(factor(Dict, m)))

    # Find the m-th root of unity over the finite ring Kpr.
    #TODO Make it deterministic instead of random.
    ξ = rand(Kpr)
    while true
        if ξ^order == Kpr(1)
            ζ = ξ^(order ÷ m)
            if ζ ≠ Kpr(1) && all(ζ .^ test .≠ Kpr(1))
                break
            end
        end
        ξ = rand(Kpr)
    end
    ζ = ξ^(order ÷ m)

    # Determine the structure of the hypercube, 
    # and compute all the factors with respect to the hypercube structure.
    dims, gens = find_generators_mod_m(m)
    facts = Vector{Vector{Int64}}(undef, k)
    local cubegen::Vector{Int64}
    local factps::zzModPolyRingElem
    comp = ZZprz(1)

    if length(dims) == 1
        cube = [k]
        cubegen = [gens[1]]

        for i = 0:cube[1]-1
            tmp = powermod(cubegen[1], cube[1] - i, m)
            facti = Kprx(1)
            for j = 0:d-1
                facti *= x - ζ^(tmp * powermod(p, j, m))
            end

            facts[i+1] = Vector{Int64}(undef, d + 1)
            factpsi = ZZprz(0)
            for j = reverse(0:d)
                facts[i+1][j+1] = coeff(coeff(facti, j).data, 0).data
                factpsi = factpsi * z + facts[i+1][j+1]
            end

            if i == 0
                factps = factpsi
            else
                comp *= factpsi
            end
        end
    elseif length(dims) == 2
        mat = zeros(Int64, 3, 2)
        mat[1, 1], mat[2, 2] = dims[1], dims[2]
        for i1 = 0:dims[1]-1, i2 = 0:dims[2]-1
            idx = (powermod(gens[1], i1, m) * powermod(gens[2], i2, m)) % m
            if idx == p
                mat[3, 1], mat[3, 2] = i1, i2
                break
            end
        end

        smat = snf(mat)
        cube = Vector{Int64}(undef, 2)
        cube[1], cube[2] = smat.S[1, 1], smat.S[2, 2]

        genmat = round.(Int64, inv(smat.V))
        cubegen = Vector{Int64}(undef, 2)
        cubegen[1] = (powermod(gens[1], genmat[1, 1], m) * powermod(gens[2], genmat[1, 2], m)) % m
        cubegen[2] = (powermod(gens[1], genmat[2, 1], m) * powermod(gens[2], genmat[2, 2], m)) % m

        for i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2
            tmp = powermod(cubegen[1], cube[1] - i1, m) * powermod(cubegen[2], cube[2] - i2, m) % m
            facti = Kprx(1)
            for j = 0:d-1
                facti *= x - ζ^(tmp * powermod(p, j, m))
            end

            facts[idx] = Vector{Int64}(undef, d + 1)
            factpsi = ZZprz(0)
            for j = reverse(0:d)
                facts[idx][j+1] = coeff(coeff(facti, j).data, 0).data
                factpsi = factpsi * z + facts[idx][j+1]
            end

            if idx == 1
                factps = factpsi
            else
                comp *= factpsi
            end
        end
    elseif length(dims) == 3
        mat = zeros(Int64, 4, 3)
        mat[1, 1], mat[2, 2], mat[3, 3] = dims[1], dims[2], dims[3]
        for i1 = 0:dims[1]-1, i2 = 0:dims[2]-1, i3 = 0:dims[3]-1
            idx = (powermod(gens[1], i1, m) * powermod(gens[2], i2, m) * powermod(gens[3], i3, m)) % m
            if idx == p
                mat[4, 1], mat[4, 2], mat[4, 3] = i1, i2, i3
                break
            end
        end

        smat = snf(mat)
        cube = Vector{Int64}(undef, 3)
        cube[1], cube[2], cube[3] = smat.S[1, 1], smat.S[2, 2], smat.S[3, 3]

        genmat = Int64.(inv(smat.V))
        cubegen = Vector{Int64}(undef, 3)
        cubegen[1] = (powermod(gens[1], genmat[1, 1], m) * powermod(gens[2], genmat[1, 2], m) * powermod(gens[3], genmat[1, 3], m)) % m
        cubegen[2] = (powermod(gens[1], genmat[2, 1], m) * powermod(gens[2], genmat[2, 2], m) * powermod(gens[3], genmat[2, 3], m)) % m
        cubegen[3] = (powermod(gens[1], genmat[3, 1], m) * powermod(gens[2], genmat[3, 2], m) * powermod(gens[3], genmat[3, 3], m)) % m

        for i3 = 0:cube[3]-1, i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2 + cube[1] * cube[2] * i3
            tmp = powermod(cubegen[1], cube[1] - i1, m) * powermod(cubegen[2], cube[2] - i2, m) * powermod(cubegen[3], cube[3] - i3, m) % m
            facti = Kprx(1)
            for j = 0:d-1
                facti *= x - ζ^(tmp * powermod(p, j, m))
            end

            facts[idx] = Vector{Int64}(undef, d + 1)
            factpsi = ZZprz(0)
            for j = reverse(0:d)
                facts[idx][j+1] = coeff(coeff(facti, j).data, 0).data
                factpsi = factpsi * z + facts[idx][j+1]
            end

            if idx == 1
                factps = factpsi
            else
                comp *= factpsi
            end
        end
    elseif length(dims) == 4
        mat = zeros(Int64, 5, 4)
        mat[1, 1], mat[2, 2], mat[3, 3], mat[4, 4] = dims[1], dims[2], dims[3], dims[4]
        for i1 = 0:dims[1]-1, i2 = 0:dims[2]-1, i3 = 0:dims[3]-1, i4 = 0:dims[4]-1
            idx = (powermod(gens[1], i1, m) * powermod(gens[2], i2, m) * powermod(gens[3], i3, m) * powermod(gens[4], i4, m)) % m
            if idx == p
                mat[5, 1], mat[5, 2], mat[5, 3], mat[5, 4] = i1, i2, i3, i4
                break
            end
        end

        smat = snf(mat)
        cube = Vector{Int64}(undef, 4)
        cube[1], cube[2], cube[3], cube[4] = smat.S[1, 1], smat.S[2, 2], smat.S[3, 3], smat.S[4, 4]

        genmat = Int64.(inv(smat.V))
        cubegen = Vector{Int64}(undef, 4)
        cubegen[1] = (powermod(gens[1], genmat[1, 1], m) * powermod(gens[2], genmat[1, 2], m) * powermod(gens[3], genmat[1, 3], m) * powermod(gens[4], genmat[1, 4], m)) % m
        cubegen[2] = (powermod(gens[1], genmat[2, 1], m) * powermod(gens[2], genmat[2, 2], m) * powermod(gens[3], genmat[2, 3], m) * powermod(gens[4], genmat[2, 4], m)) % m
        cubegen[3] = (powermod(gens[1], genmat[3, 1], m) * powermod(gens[2], genmat[3, 2], m) * powermod(gens[3], genmat[3, 3], m) * powermod(gens[4], genmat[3, 4], m)) % m
        cubegen[4] = (powermod(gens[1], genmat[4, 1], m) * powermod(gens[2], genmat[4, 2], m) * powermod(gens[3], genmat[4, 3], m) * powermod(gens[4], genmat[4, 4], m)) % m

        for i4 = 0:cube[4]-1, i3 = 0:cube[3]-1, i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2 + cube[1] * cube[2] * i3 + cube[1] * cube[2] * cube[3] * i4
            tmp = powermod(cubegen[1], cube[1] - i1, m) * powermod(cubegen[2], cube[2] - i2, m) * powermod(cubegen[3], cube[3] - i3, m) * powermod(cubegen[4], cube[4] - i4, m) % m
            facti = Kprx(1)
            for j = 0:d-1
                facti *= x - ζ^(tmp * powermod(p, j, m))
            end

            facts[idx] = Vector{Int64}(undef, d + 1)
            factpsi = ZZprz(0)
            for j = reverse(0:d)
                facts[idx][j+1] = coeff(coeff(facti, j).data, 0).data
                factpsi = factpsi * z + facts[idx][j+1]
            end

            if idx == 1
                factps = factpsi
            else
                comp *= factpsi
            end
        end
    end

    # Find the resolution of unity for the first slot.
    resol = Vector{UInt64}(undef, N)
    R, _ = residue_ring(ZZprz, factps)

    resolp = (R(comp)^(order - 1)).data * comp
    @inbounds @simd for i = 0:N-1
        resol[i+1] = coeff(resolp, i).data
    end

    cube, cubegen, facts, resol
end

find_and_save_factors(m::Int64, p::Int64, r::Int64)::Tuple{Vector{Int64},Vector{Int64},Vector{Vector{Int64}},Vector{UInt64}} = begin
    while (r + 1) * log2(p) < 62
        r += 1
    end

    cube, cubegen, factors, resol = factor_cyclotomic(m, p, r)
    save_factors("m$(m)p$(p)", cube, cubegen, factors, resol)
    cube, cubegen, factors, resol
end

save_factors(paramname::String, cube::Vector{Int64}, cubegen::Vector{Int64}, factors::Vector{Vector{Int64}}, resol::Vector{UInt64})::Nothing = begin
    path = String(@__DIR__) * "/factors.jl"
    chmod(path, 0o777)
    open(path, "a") do file
        println(file, paramname, " = (", cube, ", ", cubegen, ", ", factors, ", ", Int64.(resol), ")")
    end

    return nothing
end

load_factors(m::Int64, p::Int64, r::Int64)::Tuple{Vector{Int64},Vector{Int64},Vector{Vector{Int64}},Vector{UInt64}} = begin
    factors = open(String(@__DIR__) * "/factors.jl") do file
        for str in eachline(file)
            if occursin("m$(m)p$(p)", str)
                idx1 = last(findfirst("(", str))
                idx2 = first(findfirst("], [", str))
                idx3 = first(findfirst("[[", str))
                idx4 = last(findfirst("]]", str))
                idx5 = first(findfirst(")", str))

                cube = parse.(Int64, split(str[idx1+2:idx2-1], ","))
                cubegen = parse.(Int64, split(str[idx2+4:idx3-4], ", "))
                factors = [parse.(Int64, split(xi, ",")) for xi = split(str[idx3+2:idx4-2], "], [")]
                resol = parse.(UInt64, split(str[idx4+4:idx5-2], ","))

                return cube, cubegen, factors, resol
            end
        end
    end

    if isnothing(factors)
        @info "No factorisation for given parameters found. Factorisation will be found and saved in factors.jl."
        find_and_save_factors(m, p, r)
    else
        factors
    end
end

struct IntPackerArb <: IntPacker
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    pr::Modulus
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    resol::Vector{UInt64}
    buff::Vector{UInt64}
    cyclo_rdtor::ReductorCycloWord
    slot_rdtor::Vector{ReductorArbWord}

    function IntPackerArb(m::Integer, p::Integer, r::Integer)::IntPackerArb
        d, N = ord(p, m), totient(m)
        k = N ÷ d
        pr = Modulus(p^r)

        cube, cubegen, factors, resol = load_factors(m, p, r)
        Bred!(resol, pr)
        buff = Vector{UInt64}(undef, m)
        cyclo_rdtor = ReductorCycloWord(m, pr)

        slot_rdtor = Vector{ReductorArbWord}(undef, k)
        for i = 1:k
            slot_rdtor[i] = ReductorArbWord(N, factors[i], pr)
        end

        new(m, N, d, k, pr, cube, cubegen, resol, buff, cyclo_rdtor, slot_rdtor)
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, msg::AbstractVector{UInt64}, packer::IntPackerArb)::Nothing
    m, N, pr, cube, cubegen, resol, buff, cyclo_rdtor = packer.m, packer.N, packer.pr, packer.cube, packer.cubegen, packer.resol, packer.buff, packer.cyclo_rdtor
    
    if length(res) ≠ N
        throw(DimensionMismatch("The length of the output message should be equal to the polynomial degree."))
    end
    if length(msg) ≠ packer.k
        throw(DimensionMismatch("The length of the input message should be equal to the number of slots."))
    end

    @. buff = 0
    if length(cubegen) == 1
        for i = 0:cube[1]-1
            tmp = powermod(cubegen[1], i, m)
            buff[1] = Bred(widemul(msg[i+1], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = Bred(widemul(msg[i+1], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    elseif length(cubegen) == 2
        for i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2
            tmp = powermod(cubegen[1], i1, m) * powermod(cubegen[2], i2, m) % m
            buff[1] = Bred(widemul(msg[idx], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = Bred(widemul(msg[idx], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    elseif length(cubegen) == 3
        for i3 = 0:cube[3]-1, i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2 + cube[1] * cube[2] * i3
            tmp = powermod(cubegen[1], i1, m) * powermod(cubegen[2], i2, m) * powermod(cubegen[3], i3, m) % m
            buff[1] = Bred(widemul(msg[idx], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = Bred(widemul(msg[idx], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    elseif length(cubegen) == 4
        for i4 = cube[4] - 1, i3 = 0:cube[3]-1, i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2 + cube[1] * cube[2] * i3 + cube[1] * cube[2] * cube[3] * i4
            tmp = powermod(cubegen[1], i1, m) * powermod(cubegen[2], i2, m) * powermod(cubegen[3], i3, m) * powermod(cubegen[4], i4, m) % m
            buff[1] = Bred(widemul(msg[idx], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = Bred(widemul(msg[idx], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    end

    return nothing
end

@views function unpack_to!(res::AbstractVector{UInt64}, pt::AbstractVector{UInt64}, packer::IntPackerArb)::Nothing
    buff, k, N, slot_rdtor = packer.buff, packer.k, packer.N, packer.slot_rdtor

    if length(res) ≠ k
        throw(DimensionMismatch("The length of the output message should be equal to the number of slots."))
    end
    if length(pt) ≠ N
        throw(DimensionMismatch("The length of the input message should be equal to the polynomial degree."))
    end

    for i = eachindex(slot_rdtor)
        @. buff[1:N] = pt
        reduce!(buff[1:N], slot_rdtor[i])
        res[i] = buff[1]
    end

    return nothing
end

"""
resolution(m, p, r) returns resolutions of cyclotomic polynomial of degree m over Z_pʳ.
"""
function resolution(m::Int64, p::Int64, r::Int64)::Vector{UInt64}
    d = ord(p, m)
    N = (m - 1) ÷ d

    Kp, _ = finite_field(p, d, " ")
    ZZy, y = ZZ["y"]
    ZZpr, _ = residue_ring(ZZ, p^r)
    ZZprz, z = polynomial_ring(ZZpr, "z")

    poly = ZZprz(lift(ZZy, defining_polynomial(Kp)))
    Kpr, _ = residue_ring(ZZprz, poly)
    Kprx, x = Kpr["x"]

    # Order of the multiplicative group of the multiplicative ring.
    order = big(p)^((r - 1) * d) * (big(p)^d - 1)

    # Find the m-th root of unity over the finite ring Kpr.
    #TODO Make it deterministic instead of random.
    ξ = rand(Kpr)
    while true
        if ξ^order == Kpr(1)
            ζ = ξ^(order ÷ m)
            ζ ≠ Kpr(1) && break
        end
        ξ = rand(Kpr)
    end
    ζ = ξ^(order ÷ m)

    # Compute the first factor.
    factx = Kprx(1)
    for j = 0:d-1
        factx *= x - ζ^(powermod(p, j, m))
    end

    # Embed into the integer coefficient polynomials.
    factz = ZZprz(0)
    for j = reverse(0:d)
        factz = factz * z + coeff(coeff(factx, j).data, 0).data
    end

    # Compute its compliment.
    cyclo = ZZprz(cyclotomic(m, y))
    compz = cyclo ÷ factz

    # Compute the resolution of unity.
    R, _ = residue_ring(ZZprz, factz)
    τ = (R(compz)^(order - 1)).data * compz

    # Extract the meaningful values.
    g = primitive_root_finder(m) % Int64
    resol = Vector{UInt64}(undef, N)
    pr, idx = p^r, 1
    tmp = pr - coeff(τ, 0).data
    for i = 1:N
        resol[i] = coeff(τ, idx).data + tmp
        resol[i] ≥ pr && (resol[i] -= pr)
        idx = idx * g % m
    end

    resol
end

find_and_save_resolution(m::Int64, p::Int64, r::Int64)::Vector{UInt64} = begin
    while (r + 1) * log2(p) < 62
        r += 1
    end

    resol = resolution(m, p, r)
    save_resolution("m$(m)p$(p)", resol)
    resol
end

save_resolution(paramname::String, resol::Vector{UInt64})::Nothing = begin
    path = String(@__DIR__) * "/resolutions.jl"
    chmod(path, 0o777)
    open(path, "a") do file
        println(file, paramname, " = ", Int64.(resol))
    end

    return nothing
end

function load_resolution(m::Int64, p::Int64, r::Int64)::Vector{UInt64}
    resol = open(String(@__DIR__) * "/resolutions.jl") do file
        for str in eachline(file)
            if occursin("m$(m)p$(p)", str)
                idx1 = last(findfirst("[", str))
                idx2 = first(findfirst("]", str))
                resol = parse.(UInt64, split(str[idx1+1:idx2-1], ","))
                return resol
            end
        end
    end

    if isnothing(resol)
        @info "No resolution of unity for given parameters found. Resolution will be found and saved in resolution.jl."
        find_and_save_resolution(m, p, r)
    else
        resol
    end
end

struct IntPackerSubring <: IntPacker
    m::Int64
    p::Int64
    r::Int64
    N::Int64
    k::Int64
    pr::Modulus
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    P::Vector{Modulus}
    ntterP::Vector{CyclicNTTransformer2a3b5c7d}
    beP2pr::BasisExtender
    buff::Vector{Vector{UInt64}}
    resol::Vector{Vector{UInt64}}
    invresol::Vector{Vector{UInt64}}

    @views function IntPackerSubring(m::Integer, d::Integer, p::Integer, r::Integer)::IntPackerSubring
        ordp = ord(p, m)
        N, k = (m - 1) ÷ d, (m - 1) ÷ ordp
        pr = Modulus(p^r)

        resolution = load_resolution(m, p, r)
        Bred!(resolution, pr)

        if N < k
            @inbounds for i = 1:k÷N-1
                add_to!(resolution[1:N], resolution[1:N], resolution[i*N+1:(i+1)*N], pr)
            end
            resize!(resolution, N)

            k = N
            ordp = d
        end

        cube = [k]
        cubegen = [invmod(primitive_root_finder(m), m)]

        Plen = ceil(Int64, (log2(k) + 2r * log2(p)) / 62)
        convlen = is2a3b5c7d(k) ? k : next2a3b(2k - 1)

        P = Modulus.(find_prime_cyclic(convlen, 62, Plen))
        ntterP = CyclicNTTransformer[CyclicNTTransformer(convlen, Pi) for Pi = P]
        beP2pr = BasisExtender(P, [pr])

        buff = [Vector{UInt64}(undef, convlen) for _ = 1:Plen+1]
        resol = [Vector{UInt64}(undef, convlen) for _ = 1:Plen]
        invresol = [Vector{UInt64}(undef, convlen) for _ = 1:Plen]

        # Generate resolution
        @inbounds for i = 1:Plen
            @. resol[i][1:k] = resolution
            @. resol[i][k+1:end] = 0
            ntt!(resol[i], ntterP[i])
        end

        # Generate inverse resolution
        isodd(ordp) && circshift!(resolution, k >> 1)
        Bmul_to!(resolution, UInt64(m), resolution, pr)
        add_to!(resolution, resolution, Bred(ordp, pr), pr)
        @inbounds for i = 1:Plen
            @. invresol[i][1:k] = resolution
            @. invresol[i][k+1:end] = 0
            ntt!(invresol[i], ntterP[i])
        end

        new(m, p, r, N, k, pr, cube, cubegen, P, ntterP, beP2pr, buff, resol, invresol)
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, packer::IntPackerSubring)::Nothing
    k, N, P, beP2pr = packer.k, packer.N, packer.P, packer.beP2pr
    ntterP, buff, resol = packer.ntterP, packer.buff, packer.resol

    if k % length(a) ≠ 0
        throw(DimensionMismatch("The length of the input message should divide the number of slots."))
    end
    if length(res) ≠ N
        throw(DimensionMismatch("The length of the output message should be equal to the polynomial degree."))
    end

    alen = length(a)
    @inbounds for i = eachindex(P)
        for j = 0:k÷alen-1
            buff[i][alen*j+1] = a[1]
            @. buff[i][alen*j+2:alen*(j+1)] = a[end:-1:2]
        end
        @. buff[i][k+1:end] = 0

        ntt!(buff[i], ntterP[i])
        lazy_Bmul_to!(buff[i], buff[i], resol[i], P[i])
        intt!(buff[i], ntterP[i])
        ntterP[i].N ≠ k && add_to!(buff[i][1:k-1], buff[i][1:k-1], buff[i][k+1:2k-1], P[i])
    end
    basis_extend_to!(buff[end:end], 1:k, buff[1:end-1], 1:k, beP2pr)

    @inbounds for i = 0:N÷k-1
        @. res[i*k+1:(i+1)*k] = buff[end][1:k]
    end

    return nothing
end

@views function unpack_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, packer::IntPackerSubring)::Nothing
    k, N, P, beP2pr = packer.k, packer.N, packer.P, packer.beP2pr
    ntterP, buff, invresol = packer.ntterP, packer.buff, packer.invresol

    if length(a) ≠ N
        throw(DimensionMismatch("The length of the input message should be equal to the polynomial degree."))
    end
    if k % length(res) ≠ 0
        throw(DimensionMismatch("The length of the output message should be equal to the number of slots."))
    end

    @inbounds for i = eachindex(P)
        buff[i][1] = a[1]
        @. buff[i][2:k] = a[k:-1:2]
        @. buff[i][k+1:end] = 0

        ntt!(buff[i], ntterP[i])
        lazy_Bmul_to!(buff[i], buff[i], invresol[i], P[i])
        intt!(buff[i], ntterP[i])
        ntterP[i].N ≠ k && add_to!(buff[i][1:k-1], buff[i][1:k-1], buff[i][k+1:2k-1], P[i])
    end
    basis_extend_to!(buff[end:end], 1:k, buff[1:end-1], 1:k, beP2pr)

    reslen = length(res)
    @. res = buff[end][1:reslen]

    return nothing
end

(::Type{IntPacker})(t::Integer, param::RingParam)::IntPacker = begin
    if isa(param, CyclicParam)
        throw(DomainError("Cyclic parameters are not supported."))
    end

    fac = factor(Dict, Int64(t))
    if length(fac) ≠ 1
        throw(DomainError("Only prime power plaintext moduli are supported."))
    end

    m, p, r = param.m, first(keys(fac)), first(values(fac))

    if (r + 1) * log2(p) >= 62
        throw(DomainError("The plaintext modulus is too large."))
    end

    if ispow2(m)
        IntPackerPow2(m, p, r)
    elseif (p - 1) % m == 0
        IntPackerNTT(m, p, r)
    elseif isa(param, SubringParam)
        IntPackerSubring(m, param.d, p, r)
    else
        IntPackerArb(m, p, r)
    end
end