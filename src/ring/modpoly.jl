struct _PolyEvaluatorNTT
    param::RingParam
    Q::Modulus
    ntter::NTTransformer
    autbuff::Union{Vector{UInt64},Missing}

    _PolyEvaluatorNTT(param::RingParam, Q::Modulus) = begin
        ntter = NTTransformer(param, Q)

        if typeof(param) == CyclotomicParam
            autbuff = Vector{UInt64}(undef, ispow2(param.m) ? param.N : param.m)
        elseif typeof(param) == SubringParam
            autbuff = missing
        end

        new(param, Q, ntter, autbuff)
    end
end

_Bred(x::Union{Int64,UInt64,Int128,UInt128}, eval::_PolyEvaluatorNTT) = _Bred(x, eval.Q)
_Bred(x::UInt64, eval::_PolyEvaluatorNTT) = _Bred(x, eval.Q)
_Bred!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bred!(x, eval.Q)
_Bred_to!(res::AbstractVector, x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bred_to!(res, x, eval.Q)

_neg(x::UInt64, eval::_PolyEvaluatorNTT) = _neg(x, eval.Q)
_neg!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _neg!(x, eval.Q)
_neg_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _neg_to!(res, x, eval.Q)

_add(x::UInt64, y::UInt64, eval::_PolyEvaluatorNTT) = _add(x, y, eval.Q)

_add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::_PolyEvaluatorNTT) = begin
    if isntt
        _add_to!(res, x, y, eval.Q)
    elseif typeof(eval.ntter) == SubringNTTransformer
        _sub_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = _add(res[1], y, eval.Q)
    end
end

_add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _add_to!(res, x, y, eval.Q)

_sub(x::UInt64, y::UInt64, eval::_PolyEvaluatorNTT) = _sub(x, y, eval.Q)

_sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::_PolyEvaluatorNTT) = begin
    if isntt
        _sub_to!(res, x, y, eval.Q)
    elseif typeof(eval.ntter) == SubringNTTransformer
        _add_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = _sub(res[1], y, eval.Q)
    end
end

_sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _sub_to!(res, x, y, eval.Q)

_ntt!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _ntt!(x, eval.ntter)
_ntt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _ntt_to!(res, x, eval.ntter)
_intt!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _intt!(x, eval.ntter)
_intt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _intt_to!(res, x, eval.ntter)

_mul(x::UInt64, y::UInt64, eval::_PolyEvaluatorNTT) = _Bmul(x, y, eval.Q)
_mul_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bmul_to!(res, x, y, eval.Q)
_mul_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bmul_to!(res, x, y, eval.Q)

_muladd_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bmuladd_to!(res, x, y, eval.Q)
_muladd_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bmuladd_to!(res, x, y, eval.Q)
_mulsub_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bmulsub_to!(res, x, y, eval.Q)
_mulsub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorNTT) = _Bmulsub_to!(res, x, y, eval.Q)

_automorphism!(idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::_PolyEvaluatorNTT) = begin
    if typeof(eval.ntter) == SubringNTTransformer
        _automorphism!(idx, x, isntt, eval.ntter)
    elseif typeof(eval.ntter) == CyclotomicNTTransformerPow2
        _automorphism!(idx, x, isntt, eval.autbuff, eval.Q)
    elseif typeof(eval.ntter) == CyclotomicNTTransformerArb
        _automorphism!(idx, x, isntt, eval.autbuff, eval.ntter, eval.ntter.rdtor)
    end
end

_automorphism_to!(res::AbstractVector{UInt64}, idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::_PolyEvaluatorNTT) = begin
    @. res = x
    _automorphism!(idx, res, isntt, eval)
end

struct _PolyEvaluatorWord
    param::RingParam
    Q::Modulus
    Plen::Int
    P::Vector{Modulus}
    ntterP::Vector{NTTransformer}
    rdtor::Union{ReductorCycloWord,Missing}
    beP2Q::BasisExtender
    buff::Array{Vector{UInt64},2}

    @views function _PolyEvaluatorWord(param::RingParam, Q::Modulus)
        P = Modulus.(collect(find_prime(param, 62, 3)))
        ntterP = [NTTransformer(param, Pi) for Pi = P]
        m, N = param.m, param.N

        if typeof(param) == CyclotomicParam
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
        elseif typeof(param) == SubringParam
            Plen = ceil(Int64, (log2(m) + 2log2(Q.Q)) / 62)
            rdtor = missing
            buff = [Vector{UInt64}(undef, N) for _ = 1:3, _ = 1:2]
            beP2Q = BasisExtender(P[1:Plen], [Q])
            new(param, Q, Plen, P, ntterP, rdtor, beP2Q, buff)
        end
    end

    # Reuse the transformer, and change Q dynamically.
    @views function _PolyEvaluatorWord(eval::_PolyEvaluatorWord, Q::Modulus)
        param, P, ntterP, rdtor, buff = eval.param, eval.P, eval.ntterP, eval.rdtor, eval.buff
        m, N = param.m, param.N

        if typeof(param) == CyclotomicParam && ispow2(m)
            Plen = ceil(Int64, (trailing_zeros(N) + 2log2(Q.Q)) / 62)
        else
            Plen = ceil(Int64, (log2(m) + 2log2(Q.Q)) / 62)
        end

        beP2Q = BasisExtender(P[1:Plen], [Q])
        rdtor = ismissing(rdtor) ? missing : ReductorCycloWord(rdtor, Q)
        new(param, Q, Plen, P, ntterP, rdtor, beP2Q, buff)
    end
end

_Bred(x::Union{UInt64,Int64,UInt128,Int128}, eval::_PolyEvaluatorWord) = _Bred(x, eval.Q)
_Bred!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _Bred!(x, eval.Q)
_Bred_to!(res::AbstractVector, x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _Bred_to!(res, x, eval.Q)

_neg(x::UInt64, eval::_PolyEvaluatorWord) = _neg(x, eval.Q)
_neg!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _neg!(x, eval.Q)
_neg_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _neg_to!(res, x, eval.Q)

_add(x::UInt64, y::UInt64, eval::_PolyEvaluatorWord) = _add(x, y, eval.Q)

_add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::_PolyEvaluatorWord) = begin
    if typeof(eval.ntterP[1]) == SubringNTTransformer
        _sub_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = _add(res[1], y, eval.Q)
    end
end

_add_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _add_to!(res, x, y, eval.Q)

_sub(x::UInt64, y::UInt64, eval::_PolyEvaluatorWord) = _sub(x, y, eval.Q)

_sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::UInt64, isntt::Bool, eval::_PolyEvaluatorWord) = begin
    if typeof(eval.ntterP[1]) == SubringNTTransformer
        _add_to!(res, x, y, eval.Q)
    else
        @. res = x
        res[1] = _sub(res[1], y, eval.Q)
    end
end

_sub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _sub_to!(res, x, y, eval.Q)

_ntt!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = nothing
_ntt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = copy!(res, x)
_intt!(x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = nothing
_intt_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = copy!(res, x)

_mul(x::UInt64, y::UInt64, eval::_PolyEvaluatorWord) = _Bmul(x, y, eval.Q)
_mul_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _Bmul_to!(res, x, y, eval.Q)

@views function _mul_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord)
    Plen, P, ntterP, rdtor, beP2Q, buff = eval.Plen, eval.P, eval.ntterP, eval.rdtor, eval.beP2Q, eval.buff

    if ismissing(rdtor)
        # Either power-of-two cyclotomic ring, or subring.
        @inbounds for i = 1:Plen
            @. buff[i, 1] = x
            @. buff[i, 2] = y
            _ntt!(buff[i, 1], ntterP[i])
            _ntt!(buff[i, 2], ntterP[i])
            _lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            _intt!(buff[i, 1], ntterP[i])
        end
        basis_extend!(buff[1:1, 2], buff[:, 1], beP2Q)
        @. res = buff[1, 2]
    else
        # Arbitrary cyclotomic case.
        N = eval.param.N
        @inbounds for i = 1:Plen
            @. buff[i, 1][1:N] = x
            @. buff[i, 1][N+1:end] = 0
            @. buff[i, 2][1:N] = y
            @. buff[i, 2][N+1:end] = 0
            _ntt!(buff[i, 1], ntterP[i])
            _ntt!(buff[i, 2], ntterP[i])
            _lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            _intt!(buff[i, 1], ntterP[i])
        end
        basis_extend!(buff[1:1, 2], buff[:, 1], beP2Q)
        _reduce!(buff[1, 2], rdtor)
        @. res = buff[1, 2][1:N]
    end
end

_muladd_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _Bmuladd_to!(res, x, y, eval.Q)

@views function _muladd_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord)
    Q, Plen, P, ntterP, rdtor, beP2Q, buff = eval.Q, eval.Plen, eval.P, eval.ntterP, eval.rdtor, eval.beP2Q, eval.buff

    if ismissing(rdtor)
        # Either power-of-two cyclotomic ring, or subring.
        @inbounds for i = 1:Plen
            @. buff[i, 1] = x
            @. buff[i, 2] = y
            _ntt!(buff[i, 1], ntterP[i])
            _ntt!(buff[i, 2], ntterP[i])
            _lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            _intt!(buff[i, 1], ntterP[i])
        end
        basis_extend!(buff[1:1, 2], buff[:, 1], beP2Q)
        _add_to!(res, res, buff[1, 2], Q)
    else
        # Arbitrary cyclotomic case.
        N = eval.param.N
        @inbounds for i = 1:Plen
            @. buff[i, 1][1:N] = x
            @. buff[i, 1][N+1:end] = 0
            @. buff[i, 2][1:N] = y
            @. buff[i, 2][N+1:end] = 0
            _ntt!(buff[i, 1], ntterP[i])
            _ntt!(buff[i, 2], ntterP[i])
            _lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            _intt!(buff[i, 1], ntterP[i])
        end
        basis_extend!(buff[1:1, 2], buff[:, 1], beP2Q)
        _reduce!(buff[1, 2], rdtor)
        _add_to!(res, res, buff[1, 2][1:N], Q)
    end
end

_mulsub_to!(res::AbstractVector{UInt64}, x::UInt64, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord) = _Bmulsub_to!(res, x, y, eval.Q)

@views function _mulsub_to!(res::AbstractVector{UInt64}, x::AbstractVector{UInt64}, y::AbstractVector{UInt64}, eval::_PolyEvaluatorWord)
    Q, Plen, P, ntterP, rdtor, beP2Q, buff = eval.Q, eval.Plen, eval.P, eval.ntterP, eval.rdtor, eval.beP2Q, eval.buff

    if ismissing(rdtor)
        # Either power-of-two cyclotomic ring, or subring.
        @inbounds for i = 1:Plen
            @. buff[i, 1] = x
            @. buff[i, 2] = y
            _ntt!(buff[i, 1], ntterP[i])
            _ntt!(buff[i, 2], ntterP[i])
            _lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            _intt!(buff[i, 1], ntterP[i])
        end
        basis_extend!(buff[1:1, 2], buff[:, 1], beP2Q)
        _sub_to!(res, res, buff[1, 2], Q)
    else
        # Arbitrary cyclotomic case.
        N = eval.param.N
        @inbounds for i = 1:Plen
            @. buff[i, 1][1:N] = x
            @. buff[i, 1][N+1:end] = 0
            @. buff[i, 2][1:N] = y
            @. buff[i, 2][N+1:end] = 0
            _ntt!(buff[i, 1], ntterP[i])
            _ntt!(buff[i, 2], ntterP[i])
            _lazy_Bmul_to!(buff[i, 1], buff[i, 1], buff[i, 2], P[i])
            _intt!(buff[i, 1], ntterP[i])
        end
        basis_extend!(buff[1:1, 2], buff[:, 1], beP2Q)
        _reduce!(buff[1, 2], rdtor)
        _sub_to!(res, res, buff[1, 2][1:N], Q)
    end
end

_automorphism!(idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::_PolyEvaluatorWord) = begin
    if typeof(eval.ntterP[1]) == SubringNTTransformer
        _automorphism!(idx, x, false, eval.ntterP[1])
    elseif typeof(eval.ntterP[1]) == CyclotomicNTTransformerPow2
        _automorphism!(idx, x, false, eval.buff[1, 1], eval.Q)
    else
        _automorphism!(idx, x, false, eval.buff[1, 1], eval.ntterP[1], eval.rdtor)
    end
end

_automorphism_to!(res::AbstractVector{UInt64}, idx::Int64, x::AbstractVector{UInt64}, isntt::Bool, eval::_PolyEvaluatorWord) = begin
    @. res = x
    _automorphism!(idx, res, false, eval)
end

const _PolyEvaluator = Union{_PolyEvaluatorNTT,_PolyEvaluatorWord}

(::Type{_PolyEvaluator})(param::RingParam, Q::Modulus) = check_modulus(param, Q.Q) ? _PolyEvaluatorNTT(param, Q) : _PolyEvaluatorWord(param, Q)

#=================================================================================================================#

struct PolyEvaluator <: AbstractVector{_PolyEvaluator}
    param::RingParam
    evals::Vector{_PolyEvaluator}
    moduli::Vector{Modulus}

    PolyEvaluator(param::RingParam, moduli::Moduli) = begin
        evals = Vector{_PolyEvaluator}(undef, length(moduli))
        wordidx = findfirst(x -> !check_modulus(param, x.Q), moduli)
        baseeval = isnothing(wordidx) ? nothing : _PolyEvaluator(param, moduli[wordidx])

        @inbounds for i = eachindex(moduli)
            if check_modulus(param, moduli[i].Q)
                evals[i] = _PolyEvaluatorNTT(param, moduli[i])
            else
                evals[i] = _PolyEvaluatorWord(baseeval, moduli[i])
            end
        end

        new(param, evals, collect(moduli))
    end

    PolyEvaluator(param::RingParam, evals::AbstractVector{_PolyEvaluator}, moduli::Moduli) = begin
        new(param, collect(evals), collect(moduli))
    end
end

Base.:size(eval::PolyEvaluator) = size(eval.evals)
Base.:length(eval::PolyEvaluator) = length(eval.evals)
Base.:getindex(eval::PolyEvaluator, i::Int) = getindex(eval.evals, i)
Base.:getindex(eval::PolyEvaluator, idx::AbstractRange{Int64}) = PolyEvaluator(eval.param, eval.evals[idx], eval.moduli[idx])
Base.:vcat(eval::Vararg{Union{PolyEvaluator,_PolyEvaluator},N}) where {N} = begin
    param = eval[1].param
    evals = _PolyEvaluator[]
    moduli = Modulus[]

    @inbounds for i = 1:N
        if typeof(eval[i]) == PolyEvaluator
            @assert param == eval[i].param "Input polyevaluators should have the same parameters."
            evals = vcat(evals, eval[i].evals)
            moduli = vcat(moduli, eval[i].moduli)
        else
            @assert param == eval[i].param "Input polyevaluators should have the same parameters."
            evals = vcat(evals, eval[i])
            moduli = vcat(moduli, eval[i].Q)
        end
    end

    PolyEvaluator(param, evals, moduli)
end

export PolyEvaluator

#=================================================================================================================#

"""
ModScalar is a struct to store scalar in RNS form.
"""
struct ModScalar
    vals::Vector{UInt64}

    ModScalar(vals::AbstractVector{UInt64}) =
        new(collect(vals))

    ModScalar(len::Int64) =
        new(zeros(UInt64, len))

    ModScalar(val::Union{Int64,UInt64,Int128,UInt128}, eval::PolyEvaluator) = begin
        vals = [_Bred(val, eval[i]) for i = eachindex(eval)]

        new(vals)
    end

    ModScalar(val::BigInt, eval::PolyEvaluator) = begin
        vals = Vector{UInt64}(undef, length(eval))
        @inbounds for i = eachindex(vals)
            Qi = eval.moduli[i]
            vals[i] = mod(val, Qi.Q) % UInt64
        end

        new(vals)
    end
end

Base.:length(x::ModScalar) = length(x.vals)
Base.:getindex(x::ModScalar, i::Int) = getindex(x.vals, i)
Base.:getindex(x::ModScalar, idx::AbstractRange{Int64}) = ModScalar(x.vals[idx])
Base.:similar(x::ModScalar) = ModScalar(length(x))
Base.:copy(src::ModScalar) = ModScalar(copy(src.vals))
Base.:resize!(x::ModScalar, len::Int) = begin
    if len < length(x)
        deleteat!(x.vals, len+1:length(x))
    elseif len > length(x)
        append!(x.vals, zeros(UInt64, len - length(x)))
    end
end
Base.:firstindex(x::ModScalar) = firstindex(x.vals)
Base.:lastindex(x::ModScalar) = lastindex(x.vals)

Base.:copy!(dst::ModScalar, src::ModScalar) = begin
    @assert length(dst) == length(src) "The length of the destination scalar should be same to the input scalar."
    copy!(dst.vals, src.vals)
end

initialise!(s::ModScalar) = begin
    @. s.vals = zero(UInt64)
end

neg_to!(res::ModScalar, x::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) "the length of input and output scalars should be the same."
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        res.vals[i] = _neg(x.vals[i], eval[i])
    end
end

add_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) == length(y) "the length of input and output scalars should be the same."
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        res.vals[i] = _add(x.vals[i], y.vals[i], eval[i])
    end
end

sub_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) == length(y) "The length of input and output scalars should be the same."
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        res.vals[i] = _sub(x.vals[i], y.vals[i], eval[i])
    end
end

mul_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) == length(y) "The length of input and output scalars should be the same."
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        res.vals[i] = _mul(x.vals[i], y.vals[i], eval[i])
    end
end

muladd_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) == length(y) "The length of input and output scalars should be the same."
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        res.vals[i] = _add(res.vals[i], _mul(x.vals[i], y.vals[i], eval[i]), eval[i])
    end
end

mulsub_to!(res::ModScalar, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) == length(y) "The length of input and output scalars should be the same."
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        res.vals[i] = _sub(res.vals[i], _mul(x.vals[i], y.vals[i], eval[i]), eval[i])
    end
end

reduce!(x::ModScalar, eval::PolyEvaluator) = begin
    @assert length(x) == length(eval) "The parameters of the input scalar and evaluator should match."

    @inbounds for i = eachindex(eval)
        x.vals[i] = _Bred(x.vals[i], eval[i].Q)
    end
end

# There can be more optimisations.
function to_big(x::ModScalar, eval::PolyEvaluator)
    @assert length(x) == length(eval) "The lengths of input scalar and evaluator do not match."

    res = big(0)
    Qbig = prod(eval.moduli)

    @inbounds for i = eachindex(eval)
        Qi = eval.moduli[i]

        Qtilde = Qbig ÷ Qi.Q
        Qtildeinv = invmod(UInt64(Qtilde % Qi.Q), Qi.Q)
        res += _mul(x.vals[i], Qtildeinv, eval[i]) * Qtilde
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

    ModPoly(coeffs::Vector{Vector{UInt64}}, isntt::Bool) =
        new(collect(coeffs), length(coeffs[1]), Ref(isntt))

    ModPoly(N::Int64, len::Int64; isntt::Bool=true) =
        new([zeros(UInt64, N) for _ = 1:len], N, Ref(isntt))

    ModPoly(coeff::Vector{<:Union{Int64,UInt64}}, eval::PolyEvaluator; isntt::Bool=true) = begin
        coeffs = Vector{Vector{UInt64}}(undef, length(eval))
        @inbounds for i = eachindex(eval)
            coeffs[i] = Vector{UInt64}(undef, length(coeff))
            _Bred_to!(coeffs[i], coeff, eval[i])
            isntt && _ntt!(coeffs[i], eval[i])
        end

        new(coeffs, length(coeff), Ref(isntt))
    end

    ModPoly(coeff::Vector{BigInt}, eval::PolyEvaluator; isntt::Bool=true) = begin
        coeffs = Vector{Vector{UInt64}}(undef, length(eval))
        @inbounds for j = eachindex(eval)
            coeffs[j] = Vector{UInt64}(undef, length(coeff))
            for i = eachindex(coeff)
                coeffs[j][i] = mod(coeff[i], eval[j].Q.Q) % UInt64
            end
            isntt && _ntt!(coeffs[j], eval[j])
        end

        new(coeffs, length(coeff), Ref(isntt))
    end
end

Base.:length(x::ModPoly) = length(x.coeffs)
Base.:getindex(x::ModPoly, i::Int) = x.coeffs[i]
Base.:getindex(x::ModPoly, idx::AbstractRange{Int64}) = ModPoly(x.coeffs[idx], x.isntt[])
Base.:setindex!(x::ModPoly, xi::Vector{UInt64}, i::Int64) = x.coeffs[i] = xi
Base.:similar(x::ModPoly) = ModPoly(x.N, length(x), isntt=x.isntt[])
Base.:resize!(x::ModPoly, len::Int) = begin
    if len < length(x)
        deleteat!(x.coeffs, len+1:length(x))
    elseif len > length(x)
        @inbounds for _ = length(x)+1:len
            push!(x.coeffs, zeros(UInt64, x.N))
        end
    end
end
Base.:firstindex(x::ModPoly) = firstindex(x.coeffs)
Base.:lastindex(x::ModPoly) = lastindex(x.coeffs)

Base.:copy(x::ModPoly) = ModPoly(copy(x.coeffs), x.isntt[])

Base.:copy!(dst::ModPoly, src::ModPoly) = begin
    @assert dst.N == src.N "The ring dimension of the destination polynomial should be same to the input polynomial."
    @assert length(dst) == length(src) "The length of the destination polynomial should be same to the input polynomial."
    for i = eachindex(dst.coeffs)
        copy!(dst.coeffs[i], src.coeffs[i])
    end
    dst.isntt[] = src.isntt[]
end

initialise!(p::ModPoly; isntt::Bool=true) = begin
    @inbounds for i = eachindex(p.coeffs)
        @. p.coeffs[i] = zero(UInt64)
    end
    p.isntt[] = isntt
end

neg_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator) = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    @assert length(res) == length(x) == length(eval) "The parameters of the input and output polynomials and evaluator should match."

    @inbounds for i = eachindex(eval)
        _neg_to!(res.coeffs[i], x.coeffs[i], eval[i])
    end

    res.isntt[] = x.isntt[]
end

add_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert res.N == x.N == y.N "The polynomials should have the same ring degree."
    @assert length(res) == length(x) == length(y) "the length of input and output polynomials should be the same."
    @assert x.isntt[] == y.isntt[] "Input polynomials should be both in coeff form or NTT form."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _add_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end

    res.isntt[] = x.isntt[]
end

add_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _add_to!(res.coeffs[i], x.coeffs[i], y.vals[i], x.isntt[], eval[i])
    end

    res.isntt[] = x.isntt[]
end

add_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = add_to!(res, y, x, eval)

sub_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert res.N == x.N == y.N "The polynomials should have the same ring degree."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert x.isntt[] == y.isntt[] "Input polynomials should be both in coeff form or NTT form."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _sub_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end

    res.isntt[] = x.isntt[]
end

sub_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluator) = begin
    @assert res.N == x.N "Output polynomial should have the same ring degree to the input polynomial."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _sub_to!(res.coeffs[i], x.coeffs[i], y.vals[i], x.isntt[], eval[i])
    end

    res.isntt[] = x.isntt[]
end

sub_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    neg_to!(res, y, eval)
    add_to!(res, x, res, eval)
end

mul_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert x.N == y.N == res.N "The input and output polynomials should have the same length."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert x.isntt[] && y.isntt[] "Input polynomials should be both in NTT form."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _mul_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end

    res.isntt[] = true
end

mul_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluator) = mul_to!(res, y, x, eval)

mul_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    @assert y.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _mul_to!(res.coeffs[i], x.vals[i], y.coeffs[i], eval[i])
    end

    res.isntt[] = y.isntt[]
end

muladd_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert x.N == y.N == res.N "The input and output polynomials should have the same length."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert x.isntt[] && y.isntt[] "Input polynomials should be both in NTT form."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _muladd_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end

    res.isntt[] = true
end

muladd_to!(res::ModPoly, x::ModPoly, y::ModScalar, eval::PolyEvaluator) = muladd_to!(res, y, x, eval)

muladd_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    @assert y.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _muladd_to!(res.coeffs[i], x.vals[i], y.coeffs[i], eval[i])
    end
end

muladd_to!(res::ModPoly, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        tmp = _mul(x.vals[i], y.vals[i], eval[i])
        _add_to!(res.coeffs[i], res.coeffs[i], tmp, res.isntt[], eval[i])
    end
end

mulsub_to!(res::ModPoly, x::ModPoly, y::ModPoly, eval::PolyEvaluator) = begin
    @assert x.N == y.N == res.N "The input and output polynomials should have the same length."
    @assert length(res) == length(x) == length(y) "The length of input and output polynomials should be the same."
    @assert x.isntt[] && y.isntt[] "Input polynomials should be both in NTT form."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _mulsub_to!(res.coeffs[i], x.coeffs[i], y.coeffs[i], eval[i])
    end
end

mulsub_to!(res::ModPoly, x::ModScalar, y::ModPoly, eval::PolyEvaluator) = begin
    @assert y.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert length(res) == length(y) == length(x) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _mulsub_to!(res.coeffs[i], x.vals[i], y.coeffs[i], eval[i])
    end
end

mulsub_to!(res::ModPoly, x::ModScalar, y::ModScalar, eval::PolyEvaluator) = begin
    @assert length(res) == length(y) == length(x) "The length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        tmp = _mul(x.vals[i], y.vals[i], eval[i])
        _sub_to!(res.coeffs[i], res.coeffs[i], tmp, res.isntt[], eval[i])
    end
end

ntt(x::ModPoly, eval::PolyEvaluator) = begin
    res = similar(x)
    ntt_to!(res, x, eval)
    res
end

ntt!(x::ModPoly, eval::PolyEvaluator) = ntt_to!(x, x, eval)

function ntt_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator)
    @assert !x.isntt[] "Polynomial is already in NTT form."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert length(x) == length(res) "the length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _ntt_to!(res.coeffs[i], x.coeffs[i], eval[i])
    end

    res.isntt[] = true
end

intt(x::ModPoly, eval::PolyEvaluator) = begin
    res = similar(x)
    intt_to!(res, x, eval)
    res
end

intt!(x::ModPoly, eval::PolyEvaluator) = intt_to!(x, x, eval)

function intt_to!(res::ModPoly, x::ModPoly, eval::PolyEvaluator)
    @assert x.isntt[] "Polynomial is already in coefficient form."
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert length(x) == length(res) "the length of input and output polynomials should be the same."
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _intt_to!(res.coeffs[i], x.coeffs[i], eval[i])
    end

    res.isntt[] = false
end

function automorphism(x::ModPoly, idx::Int64, eval::PolyEvaluator)
    res = deepcopy(x)
    automorphism!(res, idx, eval)
end

function automorphism_to!(res::ModPoly, x::ModPoly, idx::Int64, eval::PolyEvaluator)
    @assert x.N == res.N "Input polynomial and output polynomial should have the same length."
    @assert length(x) == length(res) "The length of input and output polynomials should be the same."

    copy!(res, x)
    automorphism!(res, idx, eval)
end

function automorphism!(x::ModPoly, idx::Int64, eval::PolyEvaluator)
    @assert length(x) == length(eval) "The parameters of the input polynomial and evaluator should match."

    @inbounds for i = eachindex(eval)
        _automorphism!(idx, x.coeffs[i], x.isntt[], eval[i])
    end
end

# There can be more optimisations.
function to_big(x::ModPoly, eval::PolyEvaluator)
    @assert length(x) == length(eval) "The lengths of input polynomial and evaluator do not match."

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
            res[j] += _mul(tmp.coeffs[i][j], Qtildeinv, eval[i]) * Qtilde
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

function uniform_random_to!(us::UniformSampler, x::ModPoly, eval::PolyEvaluator)
    @assert length(x) == length(eval)
    uniform_random_to!(us, x.coeffs, eval.moduli)
end

#=============================================================================================================#

neg(x::T, eval::PolyEvaluator) where {T<:Union{ModScalar,ModPoly}} = begin
    res = similar(x)
    neg_to!(res, x, eval)
    res
end

add(x::T, y::S, eval::PolyEvaluator) where {T,S<:Union{ModScalar,ModPoly}} = begin
    res = typeof(x) == ModPoly ? similar(x) : similar(y)
    add_to!(res, x, y, eval)
    res
end

sub(x::T, y::S, eval::PolyEvaluator) where {T,S<:Union{ModScalar,ModPoly}} = begin
    res = typeof(x) == ModPoly ? similar(x) : similar(y)
    sub_to!(res, x, y, eval)
    res
end

mul(x::T, y::S, eval::PolyEvaluator) where {T,S<:Union{ModScalar,ModPoly}} = begin
    res = typeof(x) == ModPoly ? similar(x) : similar(y)
    mul_to!(res, x, y, eval)
    res
end

export ModScalar, ModPoly, neg_to!, add_to!, sub_to!, mul_to!, muladd_to!, mulsub_to!, reduce!
export ModPoly, ntt, ntt!, intt, intt!, automorphism, automorphism_to!, automorphism!, to_big
export neg, add, sub, mul, uniform_random_to!