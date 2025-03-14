@testset "test_cyclic" begin
    m = (2 * 3 * 5 * 7)^2

    param = CyclicParam(m)
    Q = Modulus(find_prime(param, 15)[1])
    us = UniformSampler()

    ntter = NTTransformer(param, Q)

    a = uniform_random(us, ntter.N, Q)
    b = deepcopy(a)

    HIENAA.ntt!(a, ntter)
    HIENAA.intt!(a, ntter)

    @test all(a .== b)
end

@testset "test_cyclotomic" begin
    m = 3^2 * 5 * 7 * 11
    param = CyclotomicParam(m)

    Q = Modulus(find_prime(param, 10)[1])
    us = UniformSampler()

    ntter = NTTransformer(param, Q)

    a = uniform_random(us, param.N, Q)
    b = deepcopy(a)

    HIENAA.ntt!(a, ntter)
    HIENAA.intt!(a, ntter)

    @test all(a .== b)
end

@testset "test_subring" begin
    m, d = 174763, 73
    param = SubringParam(m, d)

    Q = Modulus(find_prime(param, 3)[1])
    us = UniformSampler()

    ntter = NTTransformer(param, Q)

    a = uniform_random(us, ntter.N, Q)
    b = deepcopy(a)

    HIENAA.ntt!(a, ntter)
    HIENAA.intt!(a, ntter)

    @test all(a .== b)
end