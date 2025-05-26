abstract type PolyEvaluatorWord end

struct PolyEvaluatorNTT <: PolyEvaluatorWord
    param::RingParam
    Q::Modulus
    ntter::NTTransformer
    autbuff::Union{Vector{UInt64},Missing}

    PolyEvaluatorNTT(param::CyclotomicParam, Q::Modulus)::PolyEvaluatorNTT = begin
        ntter = NTTransformer(param, Q)
        autbuff = Vector{UInt64}(undef, ispow2(param.m) ? param.N : param.m)
        new(param, Q, ntter, autbuff)
    end

    PolyEvaluatorNTT(param::SubringParam, Q::Modulus)::PolyEvaluatorNTT = begin
        ntter = NTTransformer(param, Q)
        autbuff = missing
        new(param, Q, ntter, autbuff)
    end
end

Bred(x::Union{Int64,UInt64,Int128,UInt128}, eval::PolyEvaluatorNTT)::UInt64 = Bred(x, eval.Q)
Bred(x::UInt64, eval::PolyEvaluatorNTT)::UInt64 = Bred(x, eval.Q)
Bred!(x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bred!(x, eval.Q)
    return nothing
end
Bred_to!(res::AbstractVector, x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bred_to!(res, x, eval.Q)
    return nothing
end

neg(x::UInt64, eval::PolyEvaluatorNTT)::UInt64 = res = neg(x, eval.Q)
neg!(x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    neg!(x, eval.Q)
    return nothing
end
neg_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    neg_to!(res, x, eval.Q)
    return nothing
end

add(x::UInt64, y::UInt64, eval::PolyEvaluatorNTT)::UInt64 = add(x, y, eval.Q)
add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::PolyEvaluatorNTT)::Nothing = begin
    if isntt
        add_to!(res, x, y, eval.Q)
    elseif isa(eval.param, SubringParam)
        sub_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = add(res[1], y, eval.Q)
    end

    return nothing
end
add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT) = begin
    add_to!(res, x, y, eval.Q)
    return nothing
end

sub(x::UInt64, y::UInt64, eval::PolyEvaluatorNTT)::UInt64 = sub(x, y, eval.Q)
sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::PolyEvaluatorNTT)::Nothing = begin
    if isntt
        sub_to!(res, x, y, eval.Q)
    elseif isa(eval.param, SubringParam)
        add_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = sub(res[1], y, eval.Q)
    end

    return nothing
end
sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    sub_to!(res, x, y, eval.Q)
    return nothing
end

ntt!(x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    ntt!(x, eval.ntter)
    return nothing
end
ntt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    ntt_to!(res, x, eval.ntter)
    return nothing
end
intt!(x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    intt!(x, eval.ntter)
    return nothing
end
intt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    intt_to!(res, x, eval.ntter)
    return nothing
end

mul(x::UInt64, y::UInt64, eval::PolyEvaluatorNTT)::UInt64 = Bmul(x, y, eval.Q)
mul_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bmul_to!(res, x, y, eval.Q)
    return nothing
end
mul_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bmul_to!(res, x, y, eval.Q)
    return nothing
end

muladd_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bmuladd_to!(res, x, y, eval.Q)
    return nothing
end
muladd_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bmuladd_to!(res, x, y, eval.Q)
    return nothing
end
mulsub_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bmulsub_to!(res, x, y, eval.Q)
    return nothing
end
mulsub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorNTT)::Nothing = begin
    Bmulsub_to!(res, x, y, eval.Q)
    return nothing
end

automorphism!(idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::PolyEvaluatorNTT)::Nothing = begin
    if isa(eval.param, SubringParam)
        _automorphism!(idx, x, isntt, eval.ntter)
    elseif ispow2(eval.param.m)
        _automorphism!(idx, x, isntt, eval.autbuff, eval.Q)
    else
        _automorphism!(idx, x, isntt, eval.autbuff, eval.ntter, eval.ntter.rdtor)
    end

    return nothing
end

automorphism_to!(res::AbstractVector{UInt64}, idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::PolyEvaluatorNTT)::Nothing = begin
    @. res = x
    automorphism!(idx, res, isntt, eval)
    return nothing
end

struct PolyEvaluatorArb <: PolyEvaluatorWord
    param::RingParam
    Q::Modulus
    Plen::Int
    P::Vector{Modulus}
    ntterP::Vector{NTTransformer}
    rdtor::Union{ReductorCycloWord,Missing}
    beP2Q::BasisExtender
    buff::Array{Vector{UInt64},2}

    @views function PolyEvaluatorArb(param::RingParam, Q::Modulus)::PolyEvaluatorArb
        P = Modulus.(collect(find_prime(param, 62, 3)))
        ntterP = [NTTransformer(param, Pi) for Pi = P]
        m, N = param.m, param.N

        if isa(param, CyclotomicParam)
            if ispow2(m)
                Plen = ceil(Int64, (trailing_zeros(N) + 2log2(Q.Q)) / 62)
                rdtor = missing
                buff = [Vector{UInt64}(undef, N) for _ = 1:3, _ = 1:2]
                beP2Q = BasisExtender(P[1:Plen], [Q])
                new(param, Q, Plen, P, ntterP, rdtor, beP2Q, buff)
            else
                Plen = ceil(Int64, (log2(m) + 2log2(Q.Q)) / 62)
                rdtor = ReductorCycloWord(m, Q)
                buff = [Vector{UInt64}(undef, m) for _ = 1:3, _ = 1:2]
                beP2Q = BasisExtender(P[1:Plen], [Q])
                new(param, Q, Plen, P, ntterP, rdtor, beP2Q, buff)
            end
        elseif isa(param, SubringParam)
            Plen = ceil(Int64, (log2(m) + 2log2(Q.Q)) / 62)
            rdtor = missing
            buff = [Vector{UInt64}(undef, N) for _ = 1:3, _ = 1:2]
            beP2Q = BasisExtender(P[1:Plen], [Q])
            new(param, Q, Plen, P, ntterP, rdtor, beP2Q, buff)
        end
    end

    # Reuse the transformer, and change Q dynamically.
    @views function PolyEvaluatorArb(eval::PolyEvaluatorArb, Q::Modulus)::PolyEvaluatorArb
        param, P, ntterP, rdtor, buff = eval.param, eval.P, eval.ntterP, eval.rdtor, eval.buff
        m, N = param.m, param.N

        if isa(param, CyclotomicParam) && ispow2(m)
            Plen = ceil(Int64, (trailing_zeros(N) + 2log2(Q.Q)) / 62)
        else
            Plen = ceil(Int64, (log2(m) + 2log2(Q.Q)) / 62)
        end

        beP2Q = BasisExtender(P[1:Plen], [Q])
        rdtor = ismissing(rdtor) ? missing : ReductorCycloWord(rdtor, Q)
        new(param, Q, Plen, P, ntterP, rdtor, beP2Q, buff)
    end
end

Bred(x::Union{UInt64,Int64,UInt128,Int128}, eval::PolyEvaluatorArb)::UInt64 = Bred(x, eval.Q)
Bred!(x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    Bred!(x, eval.Q)
    return nothing
end
Bred_to!(res::AbstractVector, x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    Bred_to!(res, x, eval.Q)
    return nothing
end

neg(x::UInt64, eval::PolyEvaluatorArb)::UInt64 = neg(x, eval.Q)
neg!(x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    neg!(x, eval.Q)
    return nothing
end
neg_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    neg_to!(res, x, eval.Q)
    return nothing
end

add(x::UInt64, y::UInt64, eval::PolyEvaluatorArb)::UInt64 = add(x, y, eval.Q)
add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::PolyEvaluatorArb)::Nothing = begin
    if isa(eval.param, SubringParam)
        sub_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = add(res[1], y, eval.Q)
    end

    return nothing
end
add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    add_to!(res, x, y, eval.Q)
    return nothing
end

sub(x::UInt64, y::UInt64, eval::PolyEvaluatorArb)::UInt64 = sub(x, y, eval.Q)
sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::PolyEvaluatorArb)::Nothing = begin
    if isa(eval.param, SubringParam)
        add_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = sub(res[1], y, eval.Q)
    end

    return nothing
end
sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    sub_to!(res, x, y, eval.Q)
    return nothing
end

ntt!(x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = nothing
ntt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    copy!(res, x)
    return nothing
end
intt!(x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = nothing
intt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    copy!(res, x)
    return nothing
end

mul(x::UInt64, y::UInt64, eval::PolyEvaluatorArb)::UInt64 = Bmul(x, y, eval.Q)
mul_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    Bmul_to!(res, x, y, eval.Q)
    return nothing
end
@views function mul_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing
    Plen, P, ntterP, rdtor, beP2Q, buff = eval.Plen, eval.P, eval.ntterP, eval.rdtor, eval.beP2Q, eval.buff

    if ismissing(rdtor)
        # Either power-of-two cyclotomic ring, or subring.
        @inbounds for i = 1:Plen
            @. buff[i, 1] = x
            @. buff[i, 2] = y
            ntt!(buff[i, 1], ntterP[i])
            ntt!(buff[i, 2], ntterP[i])
            lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            intt!(buff[i, 1], ntterP[i])
        end
        basis_extend_to!(buff[1:1, 2], buff[:, 1], beP2Q)
        @. res = buff[1, 2]
    else
        # Arbitrary cyclotomic case.
        N = eval.param.N
        @inbounds for i = 1:Plen
            @. buff[i, 1][1:N] = x
            @. buff[i, 1][N+1:end] = 0
            @. buff[i, 2][1:N] = y
            @. buff[i, 2][N+1:end] = 0
            ntt!(buff[i, 1], ntterP[i])
            ntt!(buff[i, 2], ntterP[i])
            lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            intt!(buff[i, 1], ntterP[i])
        end
        basis_extend_to!(buff[1:1, 2], buff[:, 1], beP2Q)
        reduce!(buff[1, 2], rdtor)
        @. res = buff[1, 2][1:N]
    end

    return nothing
end

muladd_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    Bmuladd_to!(res, x, y, eval.Q)
    return nothing
end
@views function muladd_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing
    Q, Plen, P, ntterP, rdtor, beP2Q, buff = eval.Q, eval.Plen, eval.P, eval.ntterP, eval.rdtor, eval.beP2Q, eval.buff

    if ismissing(rdtor)
        # Either power-of-two cyclotomic ring, or subring.
        @inbounds for i = 1:Plen
            @. buff[i, 1] = x
            @. buff[i, 2] = y
            ntt!(buff[i, 1], ntterP[i])
            ntt!(buff[i, 2], ntterP[i])
            lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            intt!(buff[i, 1], ntterP[i])
        end
        basis_extend_to!(buff[1:1, 2], buff[:, 1], beP2Q)
        add_to!(res, res, buff[1, 2], Q)
    else
        # Arbitrary cyclotomic case.
        N = eval.param.N
        @inbounds for i = 1:Plen
            @. buff[i, 1][1:N] = x
            @. buff[i, 1][N+1:end] = 0
            @. buff[i, 2][1:N] = y
            @. buff[i, 2][N+1:end] = 0
            ntt!(buff[i, 1], ntterP[i])
            ntt!(buff[i, 2], ntterP[i])
            lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            intt!(buff[i, 1], ntterP[i])
        end
        basis_extend_to!(buff[1:1, 2], buff[:, 1], beP2Q)
        reduce!(buff[1, 2], rdtor)
        add_to!(res, res, buff[1, 2][1:N], Q)
    end

    return nothing
end

mulsub_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing = begin
    Bmulsub_to!(res, x, y, eval.Q)
    return nothing
end
@views function mulsub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::PolyEvaluatorArb)::Nothing
    Q, Plen, P, ntterP, rdtor, beP2Q, buff = eval.Q, eval.Plen, eval.P, eval.ntterP, eval.rdtor, eval.beP2Q, eval.buff

    if ismissing(rdtor)
        # Either power-of-two cyclotomic ring, or subring.
        @inbounds for i = 1:Plen
            @. buff[i, 1] = x
            @. buff[i, 2] = y
            ntt!(buff[i, 1], ntterP[i])
            ntt!(buff[i, 2], ntterP[i])
            lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            intt!(buff[i, 1], ntterP[i])
        end
        basis_extend_to!(buff[1:1, 2], buff[:, 1], beP2Q)
        sub_to!(res, res, buff[1, 2], Q)
    else
        # Arbitrary cyclotomic case.
        N = eval.param.N
        @inbounds for i = 1:Plen
            @. buff[i, 1][1:N] = x
            @. buff[i, 1][N+1:end] = 0
            @. buff[i, 2][1:N] = y
            @. buff[i, 2][N+1:end] = 0
            ntt!(buff[i, 1], ntterP[i])
            ntt!(buff[i, 2], ntterP[i])
            lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            intt!(buff[i, 1], ntterP[i])
        end
        basis_extend_to!(buff[1:1, 2], buff[:, 1], beP2Q)
        reduce!(buff[1, 2], rdtor)
        sub_to!(res, res, buff[1, 2][1:N], Q)
    end

    return nothing
end

automorphism!(idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::PolyEvaluatorArb)::Nothing = begin
    if isa(eval.param, SubringParam)
        _automorphism!(idx, x, false, eval.ntterP[1])
    elseif ispow2(eval.param.m)
        _automorphism!(idx, x, false, eval.buff[1, 1], eval.Q)
    else
        _automorphism!(idx, x, false, eval.buff[1, 1], eval.ntterP[1], eval.rdtor)
    end

    return nothing
end

automorphism_to!(res::AbstractVector{UInt64}, idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::PolyEvaluatorArb)::Nothing = begin
    copy!(res, x)
    automorphism!(idx, res, false, eval)
    return nothing
end

(::Type{PolyEvaluatorWord})(param::RingParam, Q::Modulus)::PolyEvaluatorWord = check_modulus(param, Q.Q) ? PolyEvaluatorNTT(param, Q) : PolyEvaluatorArb(param, Q)

#=================================================================================================================#

struct PolyEvaluatorRNS <: AbstractVector{PolyEvaluatorWord}
    param::RingParam
    evals::Vector{PolyEvaluatorWord}
    moduli::Vector{Modulus}

    PolyEvaluatorRNS(param::RingParam, moduli::Moduli)::PolyEvaluatorRNS = begin
        evals = Vector{PolyEvaluatorWord}(undef, length(moduli))
        wordidx = findfirst(x -> !check_modulus(param, x.Q), moduli)
        baseeval = isnothing(wordidx) ? nothing : PolyEvaluatorWord(param, moduli[wordidx])

        @inbounds for i = eachindex(moduli)
            if check_modulus(param, moduli[i].Q)
                evals[i] = PolyEvaluatorNTT(param, moduli[i])
            else
                evals[i] = PolyEvaluatorArb(baseeval, moduli[i])
            end
        end

        new(param, evals, collect(moduli))
    end

    PolyEvaluatorRNS(param::RingParam, evals::AbstractVector{PolyEvaluatorWord}, moduli::Moduli)::PolyEvaluatorRNS = begin
        new(param, collect(evals), collect(moduli))
    end
end

Base.:size(eval::PolyEvaluatorRNS)::Tuple{Int64} = size(eval.evals)
Base.:length(eval::PolyEvaluatorRNS)::Int64 = length(eval.evals)
Base.:getindex(eval::PolyEvaluatorRNS, i::Int)::PolyEvaluatorWord = getindex(eval.evals, i)
Base.:getindex(eval::PolyEvaluatorRNS, idx::AbstractRange{Int64})::PolyEvaluatorRNS = PolyEvaluatorRNS(eval.param, eval.evals[idx], eval.moduli[idx])
Base.:vcat(eval::Vararg{Union{PolyEvaluatorWord,PolyEvaluatorRNS}})::PolyEvaluatorRNS = begin
    param = eval[1].param
    evals = PolyEvaluatorWord[]
    moduli = Modulus[]

    @inbounds for i = eachindex(eval)
        if isa(eval[i], PolyEvaluatorRNS)
            if param ≠ eval[i].param
                throw(ErrorException("Input polyevaluators should have the same parameters."))
            end
            evals = vcat(evals, eval[i].evals)
            moduli = vcat(moduli, eval[i].moduli)
        else
            if param ≠ eval[i].param
                throw(ErrorException("Input polyevaluators should have the same parameters."))
            end
            evals = vcat(evals, eval[i])
            moduli = vcat(moduli, eval[i].Q)
        end
    end

    PolyEvaluatorRNS(param, evals, moduli)
end

#=================================================================================================================#

"""
ModScalar is a struct to store scalar in RNS form.
"""
struct ModScalar
    vals::Vector{UInt64}

    ModScalar(vals::AbstractVector{UInt64})::ModScalar =
        new(collect(vals))

    ModScalar(len::Int64)::ModScalar =
        new(zeros(UInt64, len))

    ModScalar(val::Union{Int64,UInt64,Int128,UInt128}, eval::PolyEvaluatorRNS)::ModScalar = begin
        vals = [Bred(val, eval[i]) for i = eachindex(eval)]

        new(vals)
    end

    ModScalar(val::BigInt, eval::PolyEvaluatorRNS)::ModScalar = begin
        vals = Vector{UInt64}(undef, length(eval))
        @inbounds for i = eachindex(vals)
            Qi = eval.moduli[i]
            vals[i] = mod(val, Qi.Q) % UInt64
        end

        new(vals)
    end
end

Base.:length(x::ModScalar)::Int64 = length(x.vals)
Base.:getindex(x::ModScalar, i::Int)::UInt64 = getindex(x.vals, i)
Base.:getindex(x::ModScalar, idx::AbstractRange{Int64})::ModScalar = ModScalar(x.vals[idx])
Base.:setindex!(x::ModScalar, val::UInt64, i::Int)::Nothing = begin
    setindex!(x.vals, val, i)
    return nothing
end
Base.:similar(x::ModScalar)::ModScalar = ModScalar(length(x))
Base.:copy(src::ModScalar)::ModScalar = ModScalar(copy(src.vals))
Base.:resize!(x::ModScalar, len::Int)::Nothing = begin
    if len < length(x)
        deleteat!(x.vals, len+1:length(x))
    elseif len > length(x)
        append!(x.vals, zeros(UInt64, len - length(x)))
    end

    return nothing
end
Base.:firstindex(x::ModScalar)::Int64 = firstindex(x.vals)
Base.:lastindex(x::ModScalar)::Int64 = lastindex(x.vals)

Base.:copy!(dst::ModScalar, src::ModScalar)::Nothing = begin
    if length(dst) ≠ length(src)
        throw(DimensionMismatch("The length of the destination scalar should be same to the input scalar."))
    end
    copy!(dst.vals, src.vals)
    return nothing
end

initialise!(s::ModScalar)::Nothing = begin
    @. s.vals = zero(UInt64)
    return nothing
end

neg_to!(res::ModScalar, x::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if length(res) ≠ length(x)
        throw(DimensionMismatch("The length of input and output scalars should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        res.vals[i] = neg(x.vals[i], eval[i])
    end

    return nothing
end

add_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output scalars should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        res.vals[i] = add(x.vals[i], y.vals[i], eval[i])
    end

    return nothing
end

sub_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output scalars should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        res.vals[i] = sub(x.vals[i], y.vals[i], eval[i])
    end

    return nothing
end

mul_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output scalars should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        res.vals[i] = mul(x.vals[i], y.vals[i], eval[i])
    end

    return nothing
end

muladd_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output scalars should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        res.vals[i] = add(res.vals[i], mul(x.vals[i], y.vals[i], eval[i]), eval[i])
    end

    return nothing
end

mulsub_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output scalars should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        res.vals[i] = sub(res.vals[i], mul(x.vals[i], y.vals[i], eval[i]), eval[i])
    end

    return nothing
end

reduce!(x::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input scalar and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        x.vals[i] = Bred(x.vals[i], eval[i].Q)
    end

    return nothing
end

# There can be more optimisations.
function to_big(x::ModScalar, eval::PolyEvaluatorRNS)::BigInt
    if length(x) ≠ length(eval)
        throw(DomainError("The lengths of input scalar and evaluator do not match."))
    end

    res = big(0)
    Qbig = prod(eval.moduli)

    @inbounds for i = eachindex(eval)
        Qi = eval.moduli[i]

        Qtilde = Qbig ÷ Qi.Q
        Qtildeinv = invmod(UInt64(Qtilde % Qi.Q), Qi.Q)
        res += mul(x.vals[i], Qtildeinv, eval[i]) * Qtilde
    end

    res %= Qbig
    if res > (Qbig >> 1)
        res -= Qbig
    end

    res
end

#=================================================================================================================#

"""
ModPoly is a struct to store the polynomial in RNS form.
"""
struct ModPoly
    coeffs::Vector{Vector{UInt64}}
    N::Int64
    isntt::RefBool

    ModPoly(coeffs::Vector{Vector{UInt64}}, isntt::Bool)::ModPoly =
        new(collect(coeffs), length(coeffs[1]), Ref(isntt))

    ModPoly(N::Int64, len::Int64; isntt::Bool=true)::ModPoly =
        new([zeros(UInt64, N) for _ = 1:len], N, Ref(isntt))

    ModPoly(coeff::Vector{<:Union{Int64,UInt64}}, eval::PolyEvaluatorRNS; isntt::Bool=true)::ModPoly = begin
        coeffs = Vector{Vector{UInt64}}(undef, length(eval))
        @inbounds for i = eachindex(eval)
            coeffs[i] = Vector{UInt64}(undef, length(coeff))
            Bred_to!(coeffs[i], coeff, eval[i])
            isntt && ntt!(coeffs[i], eval[i])
        end

        new(coeffs, length(coeff), Ref(isntt))
    end

    ModPoly(coeff::Vector{BigInt}, eval::PolyEvaluatorRNS; isntt::Bool=true)::ModPoly = begin
        coeffs = Vector{Vector{UInt64}}(undef, length(eval))
        @inbounds for j = eachindex(eval)
            coeffs[j] = Vector{UInt64}(undef, length(coeff))
            for i = eachindex(coeff)
                coeffs[j][i] = mod(coeff[i], eval[j].Q.Q) % UInt64
            end
            isntt && ntt!(coeffs[j], eval[j])
        end

        new(coeffs, length(coeff), Ref(isntt))
    end
end

Base.:length(x::ModPoly)::Int64 = length(x.coeffs)
Base.:getindex(x::ModPoly, i::Int)::Vector{UInt64} = x.coeffs[i]
Base.:getindex(x::ModPoly, idx::AbstractRange{Int64})::ModPoly = ModPoly(x.coeffs[idx], x.isntt[])
Base.:setindex!(x::ModPoly, xi::Vector{UInt64}, i::Int64)::Nothing = begin
    x.coeffs[i] = xi
    return nothing
end
Base.:similar(x::ModPoly)::ModPoly = ModPoly(x.N, length(x), isntt=x.isntt[])
Base.:resize!(x::ModPoly, len::Int)::Nothing = begin
    if len < length(x)
        deleteat!(x.coeffs, len+1:length(x))
    elseif len > length(x)
        @inbounds for _ = length(x)+1:len
            push!(x.coeffs, zeros(UInt64, x.N))
        end
    end
    
    return nothing
end
Base.:firstindex(x::ModPoly)::Int64 = firstindex(x.coeffs)
Base.:lastindex(x::ModPoly)::Int64 = lastindex(x.coeffs)

Base.:copy(x::ModPoly)::ModPoly = ModPoly(copy(x.coeffs), x.isntt[])

Base.:copy!(dst::ModPoly, src::ModPoly)::Nothing = begin
    if dst.N ≠ src.N
        throw(DimensionMismatch("The ring dimension of the destination polynomial should be same to the input polynomial."))
    end
    if length(dst) ≠ length(src)
        throw(DimensionMismatch("The length of the destination polynomial should be same to the input polynomial."))
    end
    for i = eachindex(dst.coeffs)
        copy!(dst.coeffs[i], src.coeffs[i])
    end
    dst.isntt[] = src.isntt[]
    
    return nothing
end

initialise!(p::ModPoly; isntt::Bool=true)::Nothing = begin
    @inbounds for i = eachindex(p.coeffs)
        @. p.coeffs[i] = zero(UInt64)
    end
    p.isntt[] = isntt
    
    return nothing
end

neg_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if res.N ≠ x.N
        throw(DimensionMismatch("Output polynomial should have the same ring degree to the input polynomial."))
    end
    if length(res) ≠ length(x) ≠ length(eval)
        throw(DimensionMismatch("The parameters of the input and output polynomials and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        neg_to!(res.coeffs[i], x.coeffs[i], eval[i])
    end
    res.isntt[] = x.isntt[]

    return nothing
end

add_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(res.N == x.N == y.N)
        throw(DimensionMismatch("The polynomials should have the same ring degree."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("the length of input and output polynomials should be the same."))
    end
    if x.isntt[] ≠ y.isntt[]
        throw(DomainError("Input polynomials should be both in coeff form or NTT form."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        add_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end
    res.isntt[] = x.isntt[]

    return nothing
end

add_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if res.N ≠ x.N
        throw(DimensionMismatch("Output polynomial should have the same ring degree to the input polynomial."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        add_to!(res.coeffs[i], x.coeffs[i], y.vals[i], x.isntt[], eval[i])
    end
    res.isntt[] = x.isntt[]

    return nothing
end

add_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    add_to!(res, y, x, eval)
    return nothing
end

sub_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(res.N == x.N == y.N)
        throw(DimensionMismatch("The polynomials should have the same ring degree."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if x.isntt[] ≠ y.isntt[]
        throw(DomainError("Input polynomials should be both in coeff form or NTT form."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        sub_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end
    res.isntt[] = x.isntt[]

    return nothing
end

sub_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if res.N ≠ x.N
        throw(DimensionMismatch("Output polynomial should have the same ring degree to the input polynomial."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        sub_to!(res.coeffs[i], x.coeffs[i], y.vals[i], x.isntt[], eval[i])
    end
    res.isntt[] = x.isntt[]

    return nothing
end

sub_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    neg_to!(res, y, eval)
    add_to!(res, x, res, eval)
    return nothing
end

mul_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(x.N == y.N == res.N)
        throw(DimensionMismatch("The input and output polynomials should have the same length."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if !(x.isntt[] == y.isntt[] == true)
        throw(DomainError("Input polynomials should be both in NTT form."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        mul_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end
    res.isntt[] = true

    return nothing
end

mul_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    mul_to!(res, y, x, eval)
    return nothing
end

mul_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if y.N ≠ res.N
        throw(DimensionMismatch("Input polynomial and output polynomial should have the same length."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        mul_to!(res.coeffs[i], x.vals[i], y.coeffs[i], eval[i])
    end
    res.isntt[] = y.isntt[]

    return nothing
end

muladd_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(x.N == y.N == res.N)
        throw(DimensionMismatch("The input and output polynomials should have the same length."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if !(x.isntt[] == y.isntt[] == true)
        throw(DomainError("Input polynomials should be both in NTT form."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        muladd_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end
    res.isntt[] = true

    return nothing
end

muladd_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    muladd_to!(res, y, x, eval)
    return nothing
end

muladd_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if y.N ≠ res.N
        throw(DimensionMismatch("Input polynomial and output polynomial should have the same length."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end
    if res.isntt[] ≠ y.isntt[]
        throw(DomainError("Input and output polynomial should be in the same form."))
    end

    @inbounds for i = eachindex(eval)
        muladd_to!(res.coeffs[i], x.vals[i], y.coeffs[i], eval[i])
    end

    return nothing
end

muladd_to!(res::ModPoly, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        tmp = mul(x.vals[i], y.vals[i], eval[i])
        add_to!(res.coeffs[i], res.coeffs[i], tmp, res.isntt[], eval[i])
    end

    return nothing
end

mulsub_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(x.N == y.N == res.N)
        throw(DimensionMismatch("The input and output polynomials should have the same length."))
    end
    if !(length(res) == length(x) == length(y))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if !(x.isntt[] == y.isntt[] == true)
        throw(DomainError("Input polynomials should be both in NTT form."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        mulsub_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end

    return nothing
end

mulsub_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    if y.N ≠ res.N
        throw(DimensionMismatch("Input polynomial and output polynomial should have the same length."))
    end
    if !(length(res) == length(y) == length(x))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end
    if res.isntt[] ≠ y.isntt[]
        throw(DomainError("Input and output polynomial should be in the same form."))
    end

    @inbounds for i = eachindex(eval)
        mulsub_to!(res.coeffs[i], x.vals[i], y.coeffs[i], eval[i])
    end

    return nothing
end

mulsub_to!(res::ModPoly, x::ModScalar, y::ModScalar, eval::PolyEvaluatorRNS)::Nothing = begin
    if !(length(res) == length(y) == length(x))
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        tmp = mul(x.vals[i], y.vals[i], eval[i])
        sub_to!(res.coeffs[i], res.coeffs[i], tmp, res.isntt[], eval[i])
    end

    return nothing
end

ntt(x::ModPoly, eval::PolyEvaluatorRNS)::ModPoly = begin
    res = similar(x)
    ntt_to!(res, x, eval)
    res
end

ntt!(x::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    ntt_to!(x, x, eval)
    return nothing
end

function ntt_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluatorRNS)::Nothing
    if x.isntt[]
        throw(DomainError("Polynomial is already in NTT form."))
    end
    if x.N ≠ res.N
        throw(DimensionMismatch("Input polynomial and output polynomial should have the same length."))
    end
    if length(x) ≠ length(res)
        throw(DimensionMismatch("the length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        ntt_to!(res.coeffs[i], x.coeffs[i], eval[i])
    end
    res.isntt[] = true

    return nothing
end

intt(x::ModPoly, eval::PolyEvaluatorRNS)::ModPoly = begin
    res = similar(x)
    intt_to!(res, x, eval)
    res
end

intt!(x::ModPoly, eval::PolyEvaluatorRNS)::Nothing = begin
    intt_to!(x, x, eval)
    return nothing
end

function intt_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluatorRNS)::Nothing
    if !x.isntt[]
        throw(DomainError("Polynomial is already in coefficient form."))
    end
    if x.N ≠ res.N
        throw(DimensionMismatch("Input polynomial and output polynomial should have the same length."))
    end
    if length(x) ≠ length(res)
        throw(DimensionMismatch("the length of input and output polynomials should be the same."))
    end
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        intt_to!(res.coeffs[i], x.coeffs[i], eval[i])
    end
    res.isntt[] = false

    return nothing
end

function automorphism(x::ModPoly, idx::Int64, eval::PolyEvaluatorRNS)::ModPoly
    res = deepcopy(x)
    automorphism!(res, idx, eval)
    res
end

function automorphism_to!(res::ModPoly, x::ModPoly, idx::Int64, eval::PolyEvaluatorRNS)::Nothing
    if x.N ≠ res.N
        throw(DimensionMismatch("Input polynomial and output polynomial should have the same length."))
    end
    if length(x) ≠ length(res)
        throw(DimensionMismatch("The length of input and output polynomials should be the same."))
    end

    copy!(res, x)
    automorphism!(res, idx, eval)

    return nothing
end

function automorphism!(x::ModPoly, idx::Int64, eval::PolyEvaluatorRNS)::Nothing
    if length(x) ≠ length(eval)
        throw(DomainError("The parameters of the input polynomial and evaluator should match."))
    end

    @inbounds for i = eachindex(eval)
        automorphism!(idx, x.coeffs[i], x.isntt[], eval[i])
    end

    return nothing
end

# There can be more optimisations.
function to_big(x::ModPoly, eval::PolyEvaluatorRNS)::Vector{BigInt}
    if length(x) ≠ length(eval)
        throw(DomainError("The lengths of input polynomial and evaluator do not match."))
    end

    tmp = similar(x)
    copy!(tmp, x)
    tmp.isntt[] && intt_to!(tmp, tmp, eval)

    N = tmp.N
    res = zeros(BigInt, N)
    Qbig = prod(eval.moduli)

    @inbounds for i = eachindex(eval)
        Qi = eval.moduli[i]

        Qtilde = Qbig ÷ Qi.Q
        Qtildeinv = invmod(UInt64(Qtilde % Qi.Q), Qi.Q)
        for j = 1:N
            res[j] += mul(tmp.coeffs[i][j], Qtildeinv, eval[i]) * Qtilde
        end
    end

    @. res %= Qbig
    @inbounds for i = 1:N
        if res[i] > (Qbig >> 1)
            res[i] -= Qbig
        end
    end

    res
end

function uniform_random_to!(us::UniformSampler, x::ModPoly, eval::PolyEvaluatorRNS)::Nothing
    if length(x) ≠ length(eval)
        throw(DomainError("The lengths of input polynomial and evaluator do not match."))
    end

    uniform_random_to!(us, x.coeffs, eval.moduli)

    return nothing
end

#=============================================================================================================#

neg(x::T, eval::PolyEvaluatorRNS) where {T<:Union{ModScalar,ModPoly}} = begin
    res = similar(x)
    neg_to!(res, x, eval)
    res
end::Union{ModScalar,ModPoly}

add(x::T, y::S, eval::PolyEvaluatorRNS) where {T,S<:Union{ModScalar,ModPoly}} = begin
    res = typeof(x) == ModPoly ? similar(x) : similar(y)
    add_to!(res, x, y, eval)
    res
end::Union{ModScalar,ModPoly}

sub(x::T, y::S, eval::PolyEvaluatorRNS) where {T,S<:Union{ModScalar,ModPoly}} = begin
    res = typeof(x) == ModPoly ? similar(x) : similar(y)
    sub_to!(res, x, y, eval)
    res
end::Union{ModScalar,ModPoly}

mul(x::T, y::S, eval::PolyEvaluatorRNS) where {T,S<:Union{ModScalar,ModPoly}} = begin
    res = typeof(x) == ModPoly ? similar(x) : similar(y)
    mul_to!(res, x, y, eval)
    res
end::Union{ModScalar,ModPoly}
