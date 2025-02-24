struct IntPackerPow2
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    pr::Modulus
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    buff::Vector{UInt64}
    ntter::CyclotomicNTTransformerPow2

    function IntPackerPow2(m::Integer, p::Integer, r::Integer)
        @assert isprime(p) && isodd(p) "Only odd prime power plaintext moduli are supported."

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
            @error "Only the primes of the form 8k+1 are supported now."
        end
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, msg::AbstractVector{UInt64}, packer::IntPackerPow2)
    @assert length(res) == packer.N "The length of the output message should be equal to the polynomial degree."
    @assert length(msg) == packer.k "The length of the input message should be equal to the number of slots."

    k, d, ntter = packer.k, packer.d, packer.ntter

    @. res = 0
    tmp, mask = 1, 2k - 1
    @inbounds for i = 0:k>>1-1
        res[(tmp>>1)*d+1] = msg[k>>1-i]
        res[(k-tmp>>1-1)*d+1] = msg[k-i]
        tmp = 5tmp & mask
    end
    scramble!(res[1:d:end], 2)
    _intt!(res[1:d:end], ntter)
end

@views function unpack_to!(res::AbstractVector{UInt64}, pt::AbstractVector{UInt64}, packer::IntPackerPow2)
    @assert length(res) == packer.k "The length of the output message should be equal to the number of slots."
    @assert length(pt) == packer.N "The length of the input message should be equal to the polynomial degree."

    k, d, buff, ntter = packer.k, packer.d, packer.buff, packer.ntter

    @. buff = pt[1:d:end]
    _ntt!(buff, ntter)
    scramble!(buff, 2)
    tmp, mask = 1, 2k - 1
    @inbounds for i = 0:k>>1-1
        res[k>>1-i] = buff[(tmp>>1)+1]
        res[k-i] = buff[(k-tmp>>1-1)+1]
        tmp = 5tmp & mask
    end
end

struct IntPackerNTT
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    pr::Modulus
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    ntter::CyclotomicNTTransformerArb

    function IntPackerNTT(m::Integer, p::Integer, r::Integer)
        @assert isprime(p) "Only prime power plaintext moduli are supported."
        N = totient(m)
        pr = Modulus(p^r)
        ntter = CyclotomicNTTransformerArb(m, pr)

        new(m, N, 1, N, pr, ntter.dims, ntter.gens, ntter)
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, msg::AbstractVector{UInt64}, packer::IntPackerNTT)
    @assert length(res) == length(msg) == packer.N "The length of the input and output messages should be equal to the polynomial degree."

    @. res = msg
    _intt!(res, packer.ntter)
end

@views function unpack_to!(res::AbstractVector{UInt64}, pt::AbstractVector{UInt64}, packer::IntPackerNTT)
    @assert length(res) == length(pt) == packer.N "The length of the input and output messages should be equal to the polynomial degree."

    @. res = pt
    _ntt!(res, packer.ntter)
end

# Many thanks to Simon Pohmann for the helpful discussions.
# We determine slot structure using an algorithm from https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8691744
function factor_cyclotomic(m::Int64, p::Int64, r::Int64)
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

find_and_save_factors(m::Int64, p::Int64, r::Int64) = begin
    while (r + 1) * log2(p) ≤ 62
        r += 1
    end

    cube, cubegen, factors, resol = factor_cyclotomic(m, p, r)
    save_factors("m$(m)p$(p)", cube, cubegen, factors, resol)
    cube, cubegen, factors, resol
end

save_factors(paramname::String, cube::Vector{Int64}, cubegen::Vector{Int64}, factors::Vector{Vector{Int64}}, resol::Vector{UInt64}) = begin
    path = String(@__DIR__) * "/factors.jl"
    chmod(path, 0o777)
    open(path, "a") do file
        println(file, paramname, " = (", cube, ", ", cubegen, ", ", factors, ", ", Int64.(resol), ")")
    end
end

load_factors(m::Int64, p::Int64, r::Int64) = begin
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

struct IntPackerArb
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

    function IntPackerArb(m::Integer, p::Integer, r::Integer)
        d, N = ord(p, m), totient(m)
        k = N ÷ d
        pr = Modulus(p^r)

        cube, cubegen, factors, resol = load_factors(m, p, r)
        _Bred!(resol, pr)
        buff = Vector{UInt64}(undef, m)
        cyclo_rdtor = ReductorCycloWord(m, pr)

        slot_rdtor = Vector{ReductorArbWord}(undef, k)
        for i = 1:k
            slot_rdtor[i] = ReductorArbWord(N, factors[i], pr)
        end

        new(m, N, d, k, pr, cube, cubegen, resol, buff, cyclo_rdtor, slot_rdtor)
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, msg::AbstractVector{UInt64}, packer::IntPackerArb)
    m, N, pr, cube, cubegen, resol, buff, cyclo_rdtor = packer.m, packer.N, packer.pr, packer.cube, packer.cubegen, packer.resol, packer.buff, packer.cyclo_rdtor

    @assert length(res) == N "The length of the output message should be equal to the polynomial degree."
    @assert length(msg) == packer.k "The length of the input message should be equal to the number of slots."

    @. buff = 0
    if length(cubegen) == 1
        for i = 0:cube[1]-1
            tmp = powermod(cubegen[1], i, m)
            buff[1] = _Bred(widemul(msg[i+1], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = _Bred(widemul(msg[i+1], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        _reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    elseif length(cubegen) == 2
        for i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2
            tmp = powermod(cubegen[1], i1, m) * powermod(cubegen[2], i2, m) % m
            buff[1] = _Bred(widemul(msg[idx], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = _Bred(widemul(msg[idx], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        _reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    elseif length(cubegen) == 3
        for i3 = 0:cube[3]-1, i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2 + cube[1] * cube[2] * i3
            tmp = powermod(cubegen[1], i1, m) * powermod(cubegen[2], i2, m) * powermod(cubegen[3], i3, m) % m
            buff[1] = _Bred(widemul(msg[idx], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = _Bred(widemul(msg[idx], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        _reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    elseif length(cubegen) == 4
        for i4 = cube[4] - 1, i3 = 0:cube[3]-1, i2 = 0:cube[2]-1, i1 = 0:cube[1]-1
            idx = 1 + i1 + cube[1] * i2 + cube[1] * cube[2] * i3 + cube[1] * cube[2] * cube[3] * i4
            tmp = powermod(cubegen[1], i1, m) * powermod(cubegen[2], i2, m) * powermod(cubegen[3], i3, m) * powermod(cubegen[4], i4, m) % m
            buff[1] = _Bred(widemul(msg[idx], resol[1]) + buff[1], pr)
            for j = 1:N-1
                buff[(tmp*j)%m+1] = _Bred(widemul(msg[idx], resol[j+1]) + buff[(tmp*j)%m+1], pr)
            end
        end

        _reduce!(buff, cyclo_rdtor)
        @. res = buff[1:N]
    end
end

@views function unpack_to!(res::AbstractVector{UInt64}, pt::AbstractVector{UInt64}, packer::IntPackerArb)
    buff, k, N, slot_rdtor = packer.buff, packer.k, packer.N, packer.slot_rdtor

    @assert length(res) == k "The length of the output message should be equal to the number of slots."
    @assert length(pt) == N "The length of the input message should be equal to the polynomial degree."

    for i = eachindex(slot_rdtor)
        @. buff[1:N] = pt
        _reduce!(buff[1:N], slot_rdtor[i])
        res[i] = buff[1]
    end
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
    while (r + 1) * log2(p) ≤ 62
        r += 1
    end

    resol = resolution(m, p, r)
    save_resolution("m$(m)p$(p)", resol)
    resol
end

save_resolution(paramname::String, resol::Vector{UInt64}) = begin
    path = String(@__DIR__) * "/resolutions.jl"
    chmod(path, 0o777)
    open(path, "a") do file
        println(file, paramname, " = ", Int64.(resol))
    end
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

struct IntPackerSubring
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

    @views function IntPackerSubring(m::Integer, d::Integer, p::Integer, r::Integer)
        ordp = ord(p, m)
        N, k = (m - 1) ÷ d, (m - 1) ÷ ordp
        pr = Modulus(p^r)

        resolution = load_resolution(m, p, r)
        _Bred!(resolution, pr)

        if N < k
            @inbounds for i = 1:k÷N-1
                _add_to!(resolution[1:N], resolution[1:N], resolution[i*N+1:(i+1)*N], pr)
            end
            resize!(resolution, N)

            k = N
            ordp = d
        end

        cube = [k]
        cubegen = [invmod(primitive_root_finder(m), m)]

        Plen = ceil(Int64, (log2(k) + 2r * log2(p)) / 62)
        convlen = is2a3b5c7d(k) ? k : next2a3b5c7d(2k - 1)

        P = Modulus.(_find_prime_cyclic(convlen, 62, Plen))
        ntterP = CyclicNTTransformer[CyclicNTTransformer(convlen, Pi) for Pi = P]
        beP2pr = BasisExtender(P, [pr])

        buff = [Vector{UInt64}(undef, convlen) for _ = 1:Plen+1]
        resol = [Vector{UInt64}(undef, convlen) for _ = 1:Plen]
        invresol = [Vector{UInt64}(undef, convlen) for _ = 1:Plen]

        # Generate resolution
        @inbounds for i = 1:Plen
            @. resol[i][1:k] = resolution
            @. resol[i][k+1:end] = 0
            _ntt!(resol[i], ntterP[i])
        end

        # Generate inverse resolution
        isodd(ordp) && circshift!(resolution, k >> 1)
        _Bmul_to!(resolution, UInt64(m), resolution, pr)
        _add_to!(resolution, resolution, _Bred(ordp, pr), pr)
        @inbounds for i = 1:Plen
            @. invresol[i][1:k] = resolution
            @. invresol[i][k+1:end] = 0
            _ntt!(invresol[i], ntterP[i])
        end

        new(m, p, r, N, k, pr, cube, cubegen, P, ntterP, beP2pr, buff, resol, invresol)
    end
end

@views function pack_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, packer::IntPackerSubring)
    k, N, P, beP2pr = packer.k, packer.N, packer.P, packer.beP2pr
    ntterP, buff, resol = packer.ntterP, packer.buff, packer.resol

    @assert k % length(a) == 0 "The length of the input message should divide the number of slots."
    @assert length(res) == N "The length of the output message should be equal to the polynomial degree."

    alen = length(a)
    @inbounds for i = eachindex(P)
        for j = 0:k÷alen-1
            buff[i][alen*j+1] = a[1]
            @. buff[i][alen*j+2:alen*(j+1)] = a[end:-1:2]
        end
        @. buff[i][k+1:end] = 0

        _ntt!(buff[i], ntterP[i])
        _lazy_Bmul_to!(buff[i], buff[i], resol[i], P[i])
        _intt!(buff[i], ntterP[i])
        ntterP[i].N ≠ k && _add_to!(buff[i][1:k-1], buff[i][1:k-1], buff[i][k+1:2k-1], P[i])
    end
    basis_extend!(buff[end:end], 1:k, buff[1:end-1], 1:k, beP2pr)

    @inbounds for i = 0:N÷k-1
        @. res[i*k+1:(i+1)*k] = buff[end][1:k]
    end
end

@views function unpack_to!(res::AbstractVector{UInt64}, a::AbstractVector{UInt64}, packer::IntPackerSubring)
    k, N, P, beP2pr = packer.k, packer.N, packer.P, packer.beP2pr
    ntterP, buff, invresol = packer.ntterP, packer.buff, packer.invresol

    @assert length(a) == N "The length of the input message should be equal to the polynomial degree."
    @assert k % length(res) == 0 "The length of the output message should be equal to the number of slots."

    @inbounds for i = eachindex(P)
        buff[i][1] = a[1]
        @. buff[i][2:k] = a[k:-1:2]
        @. buff[i][k+1:end] = 0

        _ntt!(buff[i], ntterP[i])
        _lazy_Bmul_to!(buff[i], buff[i], invresol[i], P[i])
        _intt!(buff[i], ntterP[i])
        ntterP[i].N ≠ k && _add_to!(buff[i][1:k-1], buff[i][1:k-1], buff[i][k+1:2k-1], P[i])
    end
    basis_extend!(buff[end:end], 1:k, buff[1:end-1], 1:k, beP2pr)

    reslen = length(res)
    @. res = buff[end][1:reslen]
end

"""
    IntPacker(t::Integer, param::RingParam)

IntPacker is an abstract type which supports SIMD packing modulo t.
"""
const IntPacker = Union{IntPackerPow2,IntPackerNTT,IntPackerArb,IntPackerSubring}

(::Type{IntPacker})(t::Integer, param::RingParam) = begin
    @assert typeof(param) ≠ CyclicParam "Cyclic parameters are not supported."

    fac = factor(Dict, Int64(t))
    @assert length(fac) == 1 "Only prime power plaintext moduli are supported."

    m, p, r = param.m, first(keys(fac)), first(values(fac))

    @assert (r + 1) * log2(p) ≤ 62 "The plaintext modulus is too large."

    if ispow2(m)
        IntPackerPow2(m, p, r)
    elseif (p - 1) % m == 0
        IntPackerNTT(m, p, r)
    elseif typeof(param) == SubringParam
        IntPackerSubring(m, param.d, p, r)
    else
        IntPackerArb(m, p, r)
    end
end