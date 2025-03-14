struct BGV <: HECiphertext
    val::RLWE
    level::RefInt

    function BGV(val::RLWE, level::Integer)::BGV
        new(val, Ref(UInt64(level)))
    end

    function BGV(N::Int64, len::Int64, level::Integer; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))::BGV
        new(RLWE(N, len, isntt=isntt, isPQ=isPQ, auxQ=auxQ), Ref(UInt64(level)))
    end
end

Base.:copy!(dst::BGV, src::BGV)::Nothing = begin
    if length(dst.val) ≠ length(src.val)
        throw(DimensionMismatch("The level of input and output ciphertexts should match."))
    end
    copy!(dst.val, src.val)
    dst.level[] = src.level[]

    return nothing
end

Base.:similar(x::BGV)::BGV = BGV(similar(x.val), x.level[])

Base.:getindex(x::BGV, idx::AbstractRange{Int64})::BGV = BGV(getindex(x.val, idx), x.level[])

Base.:length(x::BGV)::Int64 = length(x.val)
