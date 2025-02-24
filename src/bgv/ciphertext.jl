struct BGV
    val::RLWE
    level::RefInt

    function BGV(val::RLWE, level::Integer)
        new(val, Ref(UInt64(level)))
    end

    function BGV(N::Int64, len::Int64, level::Integer; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(RLWE(N, len, isntt=isntt, isPQ=isPQ, auxQ=auxQ), Ref(UInt64(level)))
    end
end

Base.:copy!(dst::BGV, src::BGV) = begin
    @assert length(dst.val) == length(src.val) "The level of input and output ciphertexts should match."
    copy!(dst.val, src.val)
    dst.level[] = src.level[]
end

Base.:similar(x::BGV) = BGV(similar(x.val), x.level[])

Base.:getindex(x::BGV, idx::AbstractRange{Int64}) = BGV(getindex(x.val, idx), x.level[])

Base.:length(x::BGV) = length(x.val)

export BGV