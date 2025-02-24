@testset "drop_level_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    drop_level_to!(ct, ct, 1, scheme)
    decrypt_to!(pt, ct, scheme)
    res = decode(pt, scheme)

    @test all(isapprox.(msg, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "Padd_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    msg2 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct = encrypt(pt1, scheme)

    add_to!(ct, ct, pt2, scheme)
    decrypt_to!(pt1, ct, scheme)
    res = decode(pt1, scheme)

    @test all(isapprox.(msg1 + msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "Cadd_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    msg2 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)

    ct1 = encrypt(pt1, scheme)
    ct2 = encrypt(pt2, scheme)

    add_to!(ct1, ct1, ct2, scheme)
    decrypt_to!(pt1, ct1, scheme)
    res = decode(pt1, scheme)

    @test all(isapprox.(msg1 + msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "Psub_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    msg2 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct = encrypt(pt1, scheme)

    sub_to!(ct, ct, pt2, scheme)
    decrypt_to!(pt2, ct, scheme)
    res = decode(pt2, scheme)

    @test all(isapprox.(msg1 - msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))

    sub_to!(ct, pt1, ct, scheme)
    decrypt_to!(pt1, ct, scheme)
    res = decode(pt1, scheme)

    @test all(isapprox.(msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "Csub_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    msg2 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)

    ct1 = encrypt(pt1, scheme)
    ct2 = encrypt(pt2, scheme)

    sub_to!(ct1, ct1, ct2, scheme)
    decrypt_to!(pt1, ct1, scheme)
    res = decode(pt1, scheme)

    @test all(isapprox.(msg1 - msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "mul_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    msg2 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct = encrypt(pt1, scheme)

    mul_to!(ct, ct, pt2, scheme)
    decrypt_to!(pt1, ct, scheme)
    res = decode(pt1, scheme)

    @test all(isapprox.(msg1 .* msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "Cmul_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    msg2 = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)

    ct1 = encrypt(pt1, scheme)
    ct2 = encrypt(pt2, scheme)

    rlk = relin_keygen(scheme)
    set_relinkey!(rlk, scheme)

    mul_to!(ct1, ct1, ct2, scheme)
    decrypt_to!(pt1, ct1, scheme)
    res = decode(pt1, scheme)

    @test all(isapprox.(msg1 .* msg2, res, atol=sketch.ring_param.m / sketch.scaling_factor))
end

@testset "rotate_test" begin
    m, hw, logP, logQ = 1 << 8, 32, 100, 150

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg = ComplexDF64.(rand(ComplexF64, scheme.oper.packer.k))
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    idx = (1,)
    rtk = rotate_keygen(idx, scheme)
    set_rotate_key!(idx, rtk, scheme)
    rotate_to!(ct, ct, idx, scheme)

    decrypt_to!(pt, ct, scheme)
    res = decode(pt, scheme)

    @test all(isapprox.((@view msg[1:end-1]), (@view res[2:end]), atol=sketch.ring_param.m / sketch.scaling_factor)) && isapprox(msg[end], res[1], atol=sketch.ring_param.m / sketch.scaling_factor)
end