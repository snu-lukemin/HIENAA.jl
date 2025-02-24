struct BFV
    val::RLWE
    level::RefInt

    function BFV(val::RLWE, level::Integer)
        new(val, Ref(UInt64(level)))
    end

    function BFV(N::Int64, len::Int64, level::Integer; isntt::Bool=true, isPQ::Bool=false)
        new(RLWE(N, len, isntt=isntt, isPQ=isPQ), Ref(UInt64(level)))
    end
end

Base.:copy!(dst::BFV, src::BFV) = begin
    copy!(dst.val, src.val)
    dst.level[] = src.level[]
end

Base.:similar(x::BFV) = BFV(similar(x.val), x.level[])

Base.:getindex(x::BFV, idx::AbstractRange{Int64}) = BFV(getindex(x.val, idx), x.level[])

Base.:length(x::BFV) = length(x.val)

export BFV