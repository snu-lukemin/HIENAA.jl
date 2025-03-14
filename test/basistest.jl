@testset "test_bextend" begin
    m, Qlen, Plen = 16, 8, 10

    param = CyclotomicParam(m)
    N = param.N

    Q = Modulus.(find_prime(param, 30, Qlen))
    P = Modulus.(find_prime(param, 40, Plen))

    evalQ = PolyEvaluatorRNS(param, Q)
    evalP = PolyEvaluatorRNS(param, P)

    be = BasisExtender(Q, P)
    us = UniformSampler()

    a = ModPoly(N, Qlen, isntt=false)
    b = ModPoly(N, Plen, isntt=false)

    uniform_random_to!(us, a.coeffs, Q)

    basis_extend_to!(b.coeffs, a.coeffs, be)
    @test all(mod.(to_big(a, evalQ) - to_big(b, evalP), prod(P)) .== 0)
end

@testset "test_sscale" begin
    m, Qlen, Tlen = 16, 8, 3

    param = CyclotomicParam(m)
    N = param.N

    Q = Modulus.(find_prime(param, 30, Qlen))
    T = Modulus.(find_prime(param, 20, Tlen))

    evalQ = PolyEvaluatorRNS(param, Q)
    evalT = PolyEvaluatorRNS(param, T)

    ss = SimpleScaler(Q, T)
    us = UniformSampler()

    a = ModPoly(N, Qlen, isntt=false)
    b = ModPoly(N, Tlen, isntt=false)

    uniform_random_to!(us, a.coeffs, Q)

    simple_scale_to!(b.coeffs, a.coeffs, ss)
    
    @test all(round.(BigInt, to_big(a, evalQ) / prod(Q) * prod(T)) .== to_big(b, evalT))
end

@testset "test_cscale" begin
    m, Qlen, Plen, Tlen = 16, 10, 10, 2

    param = CyclotomicParam(m)
    N = param.N

    T = Modulus.(find_prime(param, 20, Tlen))
    Q = Modulus.(find_prime(param, 30, Qlen))
    P = Modulus.(find_prime(param, 40, Plen))
    PQ = vcat(P, Q)

    evalQ = PolyEvaluatorRNS(param, Q)
    evalPQ = PolyEvaluatorRNS(param, PQ)

    cs = ComplexScaler(PQ, Q, prod(T) // prod(Q))
    us = UniformSampler()

    a = ModPoly(N, Qlen + Plen, isntt=false)
    b = ModPoly(N, Qlen, isntt=false)

    uniform_random_to!(us, a.coeffs, PQ)

    complex_scale_to!(b.coeffs, a.coeffs, cs)

    @test all(mod.(round.(BigInt, to_big(a, evalPQ) * prod(T) // prod(Q)) - to_big(b, evalQ), prod(Q)) .== 0)
end

@testset "test_rescale" begin
    m, Qlen, Plen = 16, 16, 3

    param = CyclotomicParam(m)
    N = param.N

    Q = Modulus.(find_prime(param, 30, Qlen))
    P = Modulus.(find_prime(param, 40, Plen))
    PQ = vcat(P, Q)

    evalQ = PolyEvaluatorRNS(param, Q)
    evalPQ = PolyEvaluatorRNS(param, PQ)

    ss = SimpleScaler(PQ, Q)
    us = UniformSampler()

    a = ModPoly(N, Qlen + Plen, isntt=false)
    b = ModPoly(N, Qlen, isntt=false)

    uniform_random_to!(us, a.coeffs, PQ)

    simple_scale_to!(b.coeffs, a.coeffs, ss)

    @test all(mod.(round.(BigInt, to_big(a, evalPQ) * 1 // prod(P)), prod(Q)) .== mod.(to_big(b, evalQ), prod(Q)))
end