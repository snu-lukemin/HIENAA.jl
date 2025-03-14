@testset "bgv_poly_test" begin
    m, hw, logP, logQ = 127, 32, 62, 261

    ring_param = SubringParam(m, 1)
    sketch = BGVParamSketch(ring_param, logP, logQ, 2^8)
    scheme = BGVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    rlk = relin_keygen(scheme)
    set_relinkey!(rlk, scheme)

    coeffs = PolyCoeffs([0, 0, 0, 0, 0, 0, 256-12, 0, 13], scheme)    # digit extraction polynomial.

    msg = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    res = evaluate(coeffs, ct, scheme)
    decrypt_to!(pt, res, scheme)
    out = decode(pt, scheme)

    @test all((msg - out) .% 2 .== 0)
end

@testset "bfv_poly_test" begin
    m, hw, logP, logQ = 127, 32, 0, 322

    ring_param = SubringParam(m, 1)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^8)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    rlk = relin_keygen(scheme)
    set_relinkey!(rlk, scheme)

    coeffs = PolyCoeffs([0, 0, 0, 0, 0, 0, 256-12, 0, 13], scheme)    # digit extraction polynomial.

    msg = uniform_random(us, scheme.oper.packer.k, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    res = evaluate(coeffs, ct, scheme)
    decrypt_to!(pt, res, scheme)
    out = decode(pt, scheme)

    @test all((msg - out) .% 2 .== 0)
end

@testset "ckks_poly_test" begin
    m, hw, logP, logQ = 1 << 9, 32, 62, 261

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 2^30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    rlk = relin_keygen(scheme)
    set_relinkey!(rlk, scheme)

    coeffs = PolyCoeffs([0.1, -0.2, 0.3, -0.4, 0.5], scheme)

    msg = rand(us.rng, ComplexF64, scheme.oper.packer.k)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    res = evaluate(coeffs, ct, scheme)
    decrypt_to!(pt, res, scheme)
    out = decode(pt, scheme)

    evalmsg = 0.1 .- 0.2msg .+ 0.3 * msg .^ 2 - 0.4 * msg .^ 3 + 0.5 * msg .^ 4

    tol = sketch.ring_param.m / sketch.scaling_factor
    @test all(isapprox.(evalmsg, out, atol=tol))
end

@testset "bgv_matmul_test" begin
    m, hw, logP, logQ = 127, 64, 60, 60

    ring_param = SubringParam(m, 1)
    sketch = BGVParamSketch(ring_param, logP, logQ, 2^8)
    scheme = BGVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    packlen = scheme.oper.packer.k

    msg = uniform_random(us, packlen, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    matrix = rand(0:scheme.oper.ptxt_modulus.Q-1, packlen, packlen)
    M = PlainMatrix(matrix, ct.level[], scheme)
    list = get_required_key_list(M)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    res = mul(M, ct, scheme)
    decrypt_to!(pt, res, scheme)
    out = decode(pt, scheme)

    @test all(out .== (matrix * msg) .% 2^8)
end

@testset "bfv_matmul_test" begin
    m, hw, logP, logQ = 127, 64, 0, 120

    ring_param = SubringParam(m, 1)
    sketch = BFVParamSketch(ring_param, logP, logQ, 2^8)
    scheme = BFVScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    packlen = scheme.oper.packer.k

    msg = uniform_random(us, packlen, scheme.oper.ptxt_modulus)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    matrix = rand(0:scheme.oper.ptxt_modulus.Q-1, packlen, packlen)
    M = PlainMatrix(matrix, ct.level[], scheme)
    list = get_required_key_list(M)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    res = mul(M, ct, scheme)
    decrypt_to!(pt, res, scheme)
    out = decode(pt, scheme)

    @test all(out .== (matrix * msg) .% 2^8)
end

@testset "ckks_matmul_test" begin
    m, hw, logP, logQ = 1 << 4, 4, 60, 120

    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 30)
    scheme = CKKSScheme(sketch)

    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    packlen = scheme.oper.packer.k

    msg = rand(ComplexF64, packlen)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    matrix = rand(ComplexF64, packlen, packlen)
    M = PlainMatrix(matrix, ct.level[], scheme)
    list = get_required_key_list(M)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    res = mul(M, ct, scheme)
    decrypt_to!(pt, res, scheme)
    out = decode(pt, scheme)

    tol = sketch.ring_param.m * packlen / sketch.scaling_factor
    @test all(isapprox(matrix * msg, out, atol=tol))
end