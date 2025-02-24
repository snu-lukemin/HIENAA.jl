@testset "drop_level_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 100, 100

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    drop_level_to!(ct, ct, 1, scheme)
    decrypt_to!(pt, ct, scheme)
    res = decode(pt, scheme)

    @test all(msg .== res)
end

@testset "rescale_test" begin
    m, hw, logP, logQ = 127, 32, 100, 100

    ring_param = SubringParam(m, 1)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)
    
    rescale_to!(ct, ct, scheme)
    decrypt_to!(pt, ct, scheme)
    res = decode(pt, scheme)

    @test all(msg .== res)
end

@testset "Padd_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 40, 40

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    msg2 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct = encrypt(pt1, scheme)

    add_to!(ct, ct, pt2, scheme)
    decrypt_to!(pt1, ct, scheme)
    res = decode(pt1, scheme)
    
    @test all(mod.(msg1 + msg2, scheme.oper.ptxt_modulus.Q) .== res)
end

@testset "Cadd_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 40, 40

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    msg2 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct1 = encrypt(pt1, scheme)
    ct2 = encrypt(pt2, scheme)
    
    add_to!(ct1, ct1, ct2, scheme)
    decrypt_to!(pt1, ct1, scheme)
    res = decode(pt1, scheme)

    @test all(mod.(msg1 + msg2, scheme.oper.ptxt_modulus.Q) .== res)
end

@testset "Psub_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 40, 40

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    msg2 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct = encrypt(pt1, scheme)

    sub_to!(ct, ct, pt2, scheme)
    decrypt_to!(pt1, ct, scheme)
    res = decode(pt1, scheme)

    @test all(mod.(msg1 - msg2, scheme.oper.ptxt_modulus.Q) .== res)
end

@testset "Csub_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 40, 40

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    msg2 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct1 = encrypt(pt1, scheme)
    ct2 = encrypt(pt2, scheme)

    sub_to!(ct1, ct1, ct2, scheme)
    decrypt_to!(pt1, ct1, scheme)
    res = decode(pt1, scheme)
    
    @test all(mod.(msg1 - msg2, scheme.oper.ptxt_modulus.Q) .== res)
end

@testset "Pmul_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 80, 160

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg1 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    msg2 = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt1 = encode(msg1, scheme)
    pt2 = encode(msg2, scheme)
    ct = encrypt(pt1, scheme)

    mul_to!(ct, ct, pt2, scheme)
    decrypt_to!(pt1, ct, scheme)
    res = decode(pt1, scheme)
    
    @test all(mod.(msg1 .* msg2, scheme.oper.ptxt_modulus.Q) .== res)
end

@testset "Cmul_test" begin
    m, hw, logP, logQ, σ, τ = 127, 32, 30, 70, 5, 1024

    ring_param = SubringParam(m, 1)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    schemesk = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, schemesk)
    pk = public_keygen(schemesk)
    rlk = relin_keygen(schemesk)
    
    schemepk = BFVScheme(sketch)
    set_encryptor!(pk, schemepk, σ=σ, τ=τ)
    set_relinkey!(rlk, schemepk)

    msg1 = uniform_random(us, schemepk.oper.packer.k, schemepk.oper.ptxt_modulus)
    msg2 = uniform_random(us, schemepk.oper.packer.k, schemepk.oper.ptxt_modulus)
    pt1 = encode(msg1, schemepk)
    pt2 = encode(msg2, schemepk)
    ct1 = encrypt(pt1, schemepk)
    ct2 = encrypt(pt2, schemepk)

    mul_to!(ct1, ct1, ct2, schemepk)
    decrypt_to!(pt1, ct1, schemesk)
    res = decode(pt1, schemesk)
    
    @test all(mod.(msg1 .* msg2, schemepk.oper.ptxt_modulus.Q) .== res)
end

@testset "rotate_test" begin
    m, hw, logP, logQ = 17 * 31, 32, 80, 140

    ring_param = CyclotomicParam(m)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^5)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    msg = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    idx = (1, 1)
    rtk = rotate_keygen(idx, scheme)
    set_rotate_key!(idx, rtk, scheme)
    rotate_to!(ct, ct, idx, scheme)

    decrypt_to!(pt, ct, scheme)
    res = decode(pt, scheme)

    @test @views all(msg[1:2:end] .== vcat(res[4:2:end], res[2])) && all(msg[2:2:end] .== vcat(res[3:2:end], res[1]))
end