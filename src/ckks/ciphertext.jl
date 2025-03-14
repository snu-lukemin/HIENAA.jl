struct CKKS <: HECiphertext
    val::RLWE
    level::RefInt
    scale::RefBigFloat

    function CKKS(val::RLWE, level::Integer, scale::Real)::CKKS
        new(val, Ref(UInt64(level)), Ref(BigFloat(scale, precision=192)))
    end

    function CKKS(N::Int64, len::Int64, level::Integer; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0), scale::Real=1.0)::CKKS
        new(RLWE(N, len, isntt=isntt, isPQ=isPQ, auxQ=auxQ), Ref(UInt64(level)), Ref(BigFloat(scale, precision=192)))
    end
end

Base.:copy!(dst::CKKS, src::CKKS)::Nothing = begin
    if length(dst.val) ≠ length(src.val)
        throw(DimensionMismatch("The level of input and output ciphertexts should match."))
    end
    copy!(dst.val, src.val)
    dst.level[] = src.level[]
    dst.scale[] = src.scale[]

    return nothing
end

Base.:similar(x::CKKS)::CKKS = CKKS(similar(x.val), x.level[], x.scale[])

Base.:getindex(x::CKKS, idx::AbstractRange{Int64})::CKKS = CKKS(getindex(x.val, idx), x.level[], x.scale[])

Base.:length(x::CKKS)::Int64 = length(x.val)