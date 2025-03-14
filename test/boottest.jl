@testset "bgv_c2s_test" begin
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

    c2s = HIENAA.Bgv.get_c2s_matrix(scheme.oper)

    list = get_required_key_list(c2s)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    # Coeffs2Slots
    decrypt_to!(pt, ct, scheme)
    coeff = decode(pt, scheme, ispacking=false)

    res = mul(c2s, ct, scheme)
    decrypt_to!(pt, res, scheme)
    c2sout = decode(pt, scheme)
    @test @views all(c2sout .== coeff[1:packlen])
end

@testset "bgv_s2c_test" begin
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

    s2c = HIENAA.Bgv.get_s2c_matrix(scheme.oper)

    list = get_required_key_list(s2c)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    # Slots2Coeffs
    decrypt_to!(pt, ct, scheme)
    slots = decode(pt, scheme, ispacking=true)

    res = mul(s2c, ct, scheme)
    decrypt_to!(pt, res, scheme)
    s2cout = decode(pt, scheme, ispacking=false)

    @test @views all(s2cout[1:packlen] .== msg)
end

@testset "bfv_c2s_test" begin
    m, hw, logP, logQ = 127, 64, 60, 60

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

    c2s = HIENAA.Bfv.get_c2s_matrix(scheme.oper)

    list = get_required_key_list(c2s)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    # Coeffs2Slots
    decrypt_to!(pt, ct, scheme)
    coeff = decode(pt, scheme, ispacking=false)

    res = mul(c2s, ct, scheme)
    decrypt_to!(pt, res, scheme)
    c2sout = decode(pt, scheme)
    @test @views all(c2sout .== coeff[1:packlen])
end

@testset "bfv_s2c_test" begin
    m, hw, logP, logQ = 127, 64, 60, 60

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

    s2c = HIENAA.Bfv.get_s2c_matrix(scheme.oper)

    list = get_required_key_list(s2c)
    rtk = rotate_keygen(list, scheme)
    set_rotate_key!(list, rtk, scheme)

    # Slots2Coeffs
    decrypt_to!(pt, ct, scheme)
    slots = decode(pt, scheme, ispacking=true)

    res = mul(s2c, ct, scheme)
    decrypt_to!(pt, res, scheme)
    s2cout = decode(pt, scheme, ispacking=false)

    @test @views all(s2cout[1:packlen] .== msg)
end

@testset "ckks_bootstrap_dense_full" begin
    m, hw, logP, logQ = 1 << 10, 192, 61 * 4, 1532

    # Set scheme
    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 40)
    scheme = CKKSScheme(sketch)

    boot_param = HIENAA.Ckks.CKKSBootParameters(
        gap=64, packlen=scheme.oper.packer.k, P0_bits=58, Q0_bits=57, sparse_hw=32,
        c2s_radix=4, c2s_bits=58, mod_bits=60, K=16, sine_fold=3, sine_degree=63, inverse_degree=31,
        s2c_radix=3, s2c_bits=42, isthin=false
    )

    # Set encryptor
    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    # Keygen
    bootkey = bootstrap_keygen(boot_param, scheme)
    set_bootstrapper!(bootkey, scheme)

    # Set message.
    packlen = boot_param.packlen
    msg = rand(ComplexF64, packlen)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    ct_boot = bootstrap(ct, scheme)
    decrypt_to!(pt, ct_boot, scheme)
    out = decode(pt, scheme, ispacking=true)

    @test all(isapprox(msg, out, atol=1e-6))
end

@testset "ckks_bootstrap_thin_full" begin
    m, hw, logP, logQ = 1 << 10, 192, 61 * 4, 1532

    # Set scheme
    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 40)
    scheme = CKKSScheme(sketch)

    boot_param = HIENAA.Ckks.CKKSBootParameters(
        gap=64, packlen=scheme.oper.packer.k, P0_bits=58, Q0_bits=57, sparse_hw=32,
        c2s_radix=4, c2s_bits=58, mod_bits=60, K=16, sine_fold=3, sine_degree=63, inverse_degree=31,
        s2c_radix=3, s2c_bits=42, isthin=true
    )

    # Set encryptor
    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    # Keygen
    bootkey = bootstrap_keygen(boot_param, scheme)
    set_bootstrapper!(bootkey, scheme)

    # Set message.
    packlen = boot_param.packlen
    msg = ComplexF64.(rand(Float64, packlen))
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    ct_boot = bootstrap(ct, scheme)
    decrypt_to!(pt, ct_boot, scheme)
    out = decode(pt, scheme, ispacking=true)

    @test all(isapprox(msg, out, atol=1e-6))
end

@testset "ckks_bootstrap_dense_sparse" begin
    m, hw, logP, logQ = 1 << 10, 192, 61 * 4, 1532

    # Set scheme
    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 40)
    scheme = CKKSScheme(sketch)

    boot_param = HIENAA.Ckks.CKKSBootParameters(
        gap=64, packlen=scheme.oper.packer.k >> 2, P0_bits=58, Q0_bits=57, sparse_hw=32,
        c2s_radix=4, c2s_bits=58, mod_bits=60, K=16, sine_fold=3, sine_degree=63, inverse_degree=31,
        s2c_radix=3, s2c_bits=42, isthin=false
    )

    # Set encryptor
    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    # Keygen
    bootkey = bootstrap_keygen(boot_param, scheme)
    set_bootstrapper!(bootkey, scheme)

    # Set message.
    packlen = boot_param.packlen
    msg = rand(ComplexF64, packlen)
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    ct_boot = bootstrap(ct, scheme)
    decrypt_to!(pt, ct_boot, scheme)
    out = decode(pt, scheme, ispacking=true)[1:packlen]

    @test all(isapprox(msg, out, atol=1e-6))
end

@testset "ckks_bootstrap_thin_sparse" begin
    m, hw, logP, logQ = 1 << 10, 192, 61 * 4, 1532

    # Set scheme
    ring_param = CyclotomicParam(m)
    sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 40)
    scheme = CKKSScheme(sketch)

    boot_param = HIENAA.Ckks.CKKSBootParameters(
        gap=64, packlen=scheme.oper.packer.k >> 2, P0_bits=58, Q0_bits=57, sparse_hw=32,
        c2s_radix=4, c2s_bits=58, mod_bits=60, K=16, sine_fold=3, sine_degree=63, inverse_degree=31,
        s2c_radix=3, s2c_bits=42, isthin=true
    )

    # Set encryptor
    us = UniformSampler()
    key = ternary_ringkey(us, ring_param.N, hw)
    set_encryptor!(key, scheme)

    # Keygen
    bootkey = bootstrap_keygen(boot_param, scheme)
    set_bootstrapper!(bootkey, scheme)

    # Set message.
    packlen = boot_param.packlen
    msg = ComplexF64.(rand(Float64, packlen))
    pt = encode(msg, scheme)
    ct = encrypt(pt, scheme)

    ct_boot = bootstrap(ct, scheme)
    decrypt_to!(pt, ct_boot, scheme)
    out = decode(pt, scheme, ispacking=true)[1:packlen]

    @test all(isapprox(msg, out, atol=1e-6))
end