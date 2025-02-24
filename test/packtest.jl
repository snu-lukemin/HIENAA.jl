@testset "test_intpacker_arb" begin
    m, p, r = 13107, 2, 4
    param = CyclotomicParam(m)

    packer = IntPacker(p^r, param)
    eval = HIENAA._PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    HIENAA._mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

@testset "test_intpacker_ntt" begin
    m, r = 3 * 5 * 7, 2
    param = CyclotomicParam(m)

    p = find_prime(param, 1)[1]
    packer = IntPacker(p^r, param)
    eval = HIENAA._PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    HIENAA._mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

@testset "test_intpacker_subring" begin
    m, p, r = 8191, 2, 32
    param = SubringParam(m, HIENAA.ord(p, m) << 1)

    packer = IntPacker(p^r, param)
    eval = HIENAA._PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    HIENAA._mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

@testset "test_intpacker_nttpow2" begin
    m, p, r = 1 << 6, 17, 1
    param = CyclotomicParam(m)

    packer = IntPacker(p^r, param)
    eval = HIENAA._PolyEvaluatorWord(param, packer.pr)

    us = UniformSampler()

    res1 = Vector{UInt64}(undef, packer.N)
    res2 = Vector{UInt64}(undef, packer.N)
    msg1 = uniform_random(us, packer.k, packer.pr)
    msg2 = uniform_random(us, packer.k, packer.pr)

    out = similar(msg1)

    pack_to!(res1, msg1, packer)
    pack_to!(res2, msg2, packer)
    HIENAA._mul_to!(res1, res1, res2, eval)

    unpack_to!(out, res1, packer)

    @test all((msg1 .* msg2 .% packer.pr.Q) .== out)
end

@testset "test_complexpacker_pow2" begin
    m = 1 << 6
    param = CyclotomicParam(m)
    packer = ComplexPacker(param)

    res = Vector{Int128}(undef, packer.N)
    msg = rand(ComplexDF64, packer.N >> 1)
    out = similar(msg)
    Δ = 1 << 40

    pack_to!(res, msg, Δ, packer)
    unpack_to!(out, res, Δ, packer)

    @test all(isapprox.(msg, out, atol=packer.N / Δ))
end