struct CKKS
    val::RLWE
    level::RefInt

    function CKKS(val::RLWE, level::Integer)
        new(val, Ref(UInt64(level)))
    end

    function CKKS(N::Int64, len::Int64, level::Integer; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(RLWE(N, len, isntt=isntt, isPQ=isPQ, auxQ=auxQ), Ref(UInt64(level)))
    end
end

Base.:copy!(dst::CKKS, src::CKKS) = begin
    @assert length(dst.val) == length(src.val) "The level of input and output ciphertexts should match."
    copy!(dst.val, src.val)
    dst.level[] = src.level[]
end

Base.:similar(x::CKKS) = CKKS(similar(x.val), x.level[])

Base.:getindex(x::CKKS, idx::AbstractRange{Int64}) = CKKS(getindex(x.val, idx), x.level[])

Base.:length(x::CKKS) = length(x.val)

export CKKS