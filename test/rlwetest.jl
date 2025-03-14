@testset "test_scale" begin
    m, d, hw, logQ, σ = 65537, 32, 192, 62 * 4, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, 0, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level)
    scale_to!(ct, ct, level - 2, oper)
    res = phase(ct, entor)

    @test all(abs.(to_big(res, oper)) .< m * hw + 6σ)
end

@testset "test_scale_P" begin
    m, hw, logP, logQ, σ = 3^4 * 5, 16, 150, 62 * 4, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, isPQ=true, auxQ=UInt64(1 << 32))
    divide_by_P_to!(ct, ct, oper)
    res = phase(ct, entor)

    @test all(abs.(to_big(res, oper)) .< m * hw + 6σ)
end

@testset "test_gadgetprod" begin
    m, d, hw, logP, logQ, σ = 65537, 64, 32, 80, 250, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    Qlen = length(param.Q) - 2
    evalQ = geteval_at(Qlen, oper; auxQ=auxQ)

    a = PlainPoly(ring_param.N, Qlen, auxQ=auxQ)
    uniform_random_to!(us, a.val, evalQ)

    rlev = rlev_encrypt(PlainConst(ModScalar(1, evalQ)), entor)
    ct = gadgetprod(a, rlev, oper)

    res = phase(ct, entor)
    @test all(abs.(to_big(res, oper) - to_big(a, oper)) .< m * hw + 6σ)
end

@testset "test_relin" begin
    m, d, hw, logP, logQ, σ = 65537, 64, 32, 80, 250, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 2
    ct = rlwe_sample(entor, level, auxQ=auxQ)
    ct2 = rlwe_sample(entor, level, auxQ=auxQ)
    evalQ = geteval_at(level, oper; auxQ=auxQ)

    vals = [mul(ct.b, ct2.b, evalQ), add(mul(ct.b, ct2.a, evalQ), mul(ct.a, ct2.b, evalQ), evalQ), mul(ct.a, ct2.a, evalQ)]
    ct3 = Tensor(vals, auxQ=auxQ)

    rlk = relin_keygen(entor)
    relinearise_to!(ct, ct3, rlk, oper)

    res = phase(ct, entor)
    @test all(abs.(to_big(res, oper) .< m * hw + 6σ))
end

@testset "test_keyswitch" begin
    m, d, hw, logP, logQ, σ = 174763, 19, 32, 120, 250, 3.2

    ring_param = SubringParam(m, d)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    key2 = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)
    entor2 = SKEncryptor(key2, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, auxQ=auxQ)

    Plen = length(param.P)
    keypt = PlainPoly(entor.keyPQ[Plen+1:end])
    key2pt = PlainPoly(entor2.keyPQ[Plen+1:end])
    ksk = rlev_encrypt(keypt, entor2)
    ksk2 = rlev_encrypt(key2pt, entor)

    keyswitch_to!(ct, ct, ksk, oper)
    adec = decompose(PlainPoly(ct.a, auxQ=auxQ), oper)
    ct2 = hoisted_keyswitch(adec, ct, ksk2, oper)
    
    res1 = phase(ct, entor2)
    res2 = phase(ct2, entor)

    @test all(abs.(to_big(res1, oper)) .< m * hw + 6σ) && all(abs.(to_big(res2, oper)) .< m * hw + 6σ)
end

@testset "test_automorphism" begin
    m, hw, logP, logQ, σ = 3 * 5 * 7 * 11, 32, 120, 250, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 1
    ct = rlwe_sample(entor, level, auxQ=auxQ)

    idx = 4
    atk = automorphism_keygen(idx, entor)
    automorphism_to!(ct, ct, idx, atk, oper)
    res = phase(ct, entor)
    @test all(abs.(to_big(res, oper) .< m * hw + 6σ))
end

@testset "test_external" begin
    m, hw, logP, logQ, σ = 127 * 17, 32, 120, 250, 3.2

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    entor = SKEncryptor(key, σ, oper)

    auxQ = UInt64(1 << 32)
    level = length(param.Q) - 1
    evalQ = geteval_at(level, oper, auxQ=auxQ)

    message = ModScalar(prod(evalQ.moduli) >> 4, evalQ)
    mu = ModScalar(2, evalQ)

    ct = rlwe_encrypt(PlainConst(message, auxQ=auxQ), entor)
    rgsw = rgsw_encrypt(PlainConst(mu), entor)
    extprod_to!(ct, ct, rgsw, oper)
    res = to_big(phase(ct, entor), oper)

    @test @views round(Int64, res[1] / prod(evalQ.moduli) * 16) == 2 && all(round.(Int64, res[2:end] / prod(evalQ.moduli) * 16) .== 0)
end