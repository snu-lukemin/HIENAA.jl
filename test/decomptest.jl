@testset "test_decomp" begin
    m, logQ, logP, dlen = 16, 300, 80, 3

    ring_param = CyclotomicParam(m)
    sketch = RLWEParamSketch(ring_param, logP, logQ, dlen)
    param = RLWEParameters(sketch)
    oper = Operator(param)

    evalQ, evalP = oper.evalQ[1:end-1], oper.evalP
    evalPQ = vcat(evalP, evalQ)

    decer = oper.decer
    us = UniformSampler()

    N, Plen, Qlen = ring_param.N, length(evalP), length(evalQ)
    a = PlainPoly(N, Qlen, isntt=false)
    uniform_random_to!(us, a.val, evalQ)

    adec = HIENAA._decompose(a, decer)
    res = PlainPoly(N, Plen + Qlen, isntt=false, isPQ=true)

    len, _ = size(adec)
    for i = 1:len
        muladd_to!(res.val, decer.gvec[i][1:Plen+Qlen], adec[i], evalPQ)
    end

    @test all(to_big(a, oper) .== to_big(res, oper) .÷ prod(evalP.moduli))
end