"""
    PlainConst is a struct for plaintext constant.
"""
struct PlainConst
    val::ModScalar
    isPQ::RefBool
    auxQ::RefInt

    function PlainConst(val::ModScalar; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(val, Ref(isPQ), Ref(auxQ))
    end

    function PlainConst(len::Int64; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(ModScalar(len), Ref(isPQ), Ref(auxQ))
    end
end

Base.:copy(x::PlainConst) = PlainConst(copy(x.val), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:copy!(dst::PlainConst, src::PlainConst) = begin
    copy!(dst.val, src.val)
    dst.isPQ[] = src.isPQ[]
    dst.auxQ[] = src.auxQ[]
end
Base.:length(x::PlainConst) = length(x.val)
Base.:getindex(x::PlainConst, idx::AbstractRange{Int64}) = PlainConst(x.val[idx], isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:similar(x::PlainConst) = PlainConst(similar(x.val), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:resize!(x::PlainConst, len::Int64) = resize!(x.val, len)

initialise!(x::PlainConst; isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    initialise!(x.val)

    x.isPQ[] = isPQ
    x.auxQ[] = auxQ
end

"""
    PlainPoly is a struct for plaintext polynomial.
"""
struct PlainPoly
    val::ModPoly
    isPQ::RefBool
    auxQ::RefInt

    function PlainPoly(val::ModPoly; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(val, Ref(isPQ), Ref(auxQ))
    end

    function PlainPoly(N::Int64, len::Int64; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(ModPoly(N, len, isntt=isntt), Ref(isPQ), Ref(auxQ))
    end
end

Base.:copy(x::PlainPoly) = PlainPoly(copy(x.val), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:copy!(dst::PlainPoly, src::PlainPoly) = begin
    copy!(dst.val, src.val)
    dst.isPQ[] = src.isPQ[]
    dst.auxQ[] = src.auxQ[]
end
Base.:length(x::PlainPoly) = length(x.val)
Base.:getindex(x::PlainPoly, idx::AbstractRange{Int64}) = PlainPoly(x.val[idx], isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:similar(x::PlainPoly) = PlainPoly(similar(x.val), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:resize!(x::PlainPoly, len::Int64) = resize!(x.val, len)

initialise!(x::PlainPoly; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    initialise!(x.val, isntt=isntt)

    x.isPQ[] = isPQ
    x.auxQ[] = auxQ
end

const PlainText = Union{PlainConst,PlainPoly}

#=================================================================================================#

"""
RLWE is a struct for RLWE ciphertext.
"""
struct RLWE
    b::ModPoly
    a::ModPoly
    isPQ::RefBool
    auxQ::RefInt

    function RLWE(b::ModPoly, a::ModPoly; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        @assert b.N == a.N && length(b) == length(a) && b.isntt[] == a.isntt[] "The mask and body should have the same parameters."
        new(b, a, Ref(isPQ), Ref(auxQ))
    end

    function RLWE(b::ModPoly; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        a = similar(b)
        new(b, a, Ref(isPQ), Ref(auxQ))
    end

    function RLWE(N::Int64, len::Int64; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        new(ModPoly(N, len, isntt=isntt), ModPoly(N, len, isntt=isntt), Ref(isPQ), Ref(auxQ))
    end
end

Base.:copy(x::RLWE) = RLWE(copy(x.b), copy(x.a), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:copy!(dst::RLWE, src::RLWE) = begin
    @assert length(dst.b) == length(src.b) "The length of input and output ciphertexts should match."
    copy!(dst.b, src.b)
    copy!(dst.a, src.a)
    dst.isPQ[] = src.isPQ[]
    dst.auxQ[] = src.auxQ[]
end

Base.:length(x::RLWE) = begin
    @assert length(x.a) == length(x.b) "Something is wrong with the ciphertext."
    length(x.b)
end

Base.:getindex(x::RLWE, idx::AbstractRange{Int64}) = RLWE(x.b[idx], x.a[idx], isPQ=x.isPQ[], auxQ=length(x) ∈ idx ? x.auxQ[] : UInt64(0))
Base.:similar(x::RLWE) = RLWE(similar(x.b), similar(x.a), isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:resize!(x::RLWE, len::Int64) = begin
    resize!(x.b, len)
    resize!(x.a, len)
end

initialise!(x::RLWE; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    initialise!(x.b, isntt=isntt)
    initialise!(x.a, isntt=isntt)

    x.isPQ[] = isPQ
    x.auxQ[] = auxQ
end

#=================================================================================================#

"""
Tensor is a struct for Tensored RLWE ciphertext.
"""
struct Tensor
    vals::Vector{ModPoly}
    isPQ::RefBool
    auxQ::RefInt

    Tensor(vals::Vector{ModPoly}; isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = new(vals, Ref(isPQ), Ref(auxQ))

    function Tensor(val::ModPoly, degree::Int64; isPQ::Bool=false, auxQ::UInt64=UInt64(0))
        vals = Vector{ModPoly}(undef, degree)

        vals[1] = val
        @inbounds for i = 2:degree
            vals[i] = similar(val)
        end

        new(vals, Ref(isPQ), Ref(auxQ))
    end

    Tensor(N::Int64, len::Int64, degree::Int64=3; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = new([ModPoly(N, len, isntt=isntt) for _ = 1:degree], Ref(isPQ), Ref(auxQ))
end

Base.:axes(x::Tensor) = begin
    degree, len = length(x.vals), length(x.vals[1])
    @inbounds for i = 2:degree
        @assert len == length(x.vals[i]) "Something is wrong with the ciphertext."
    end
    (1:degree, 1:len)
end

Base.:axes(x::Tensor, i::Int64) = begin
    degree, len = size(x)
    if i == 1
        1:degree
    else
        @inbounds for j = 2:degree
            @assert len == length(x.vals[j]) "Something is wrong with the ciphertext."
        end
        1:len
    end
end

Base.:getindex(ct::Tensor, idx1::AbstractRange{Int64}, idx2::AbstractRange{Int64}) = begin
    M = length(idx1)
    @assert M ≤ length(ct.vals) "The number of indices should be less than or equal to the degree of the tensor."
    Tensor([ct.vals[i][idx2] for i = idx1], isPQ=ct.isPQ[], auxQ=ct.auxQ[])
end

Base.:getindex(ct::Tensor, idx1::Int64, idx2::AbstractRange{Int64}) = ct.vals[idx1][idx2]

Base.:getindex(ct::Tensor, idx1::Int64, idx2::Int64) = ct.vals[idx1][idx2]

Base.:getindex(ct::Tensor, idx1::Int64) = ct.vals[idx1]

Base.:size(ct::Tensor) = begin
    degree, len = length(ct.vals), length(ct.vals[1])
    @inbounds for i = 2:degree
        @assert len == length(ct.vals[i]) "Something is wrong with the ciphertext."
    end
    (degree, len)
end

Base.:copy(src::Tensor) = Tensor([copy(vali) for vali = src.vals], isPQ=src.isPQ[], auxQ=src.auxQ[])

Base.:copy!(dst::Tensor, src::Tensor) = begin
    degree, _ = size(dst)
    for i = 1:degree
        copy!(dst.vals[i], src.vals[i])
    end
    dst.isPQ[] = src.isPQ[]
    dst.auxQ[] = src.auxQ[]
end

Base.:similar(x::Tensor) = Tensor([similar(vali) for vali = x.val], isPQ=x.isPQ[], auxQ=x.auxQ[])
Base.:resize!(x::Tensor, len::Int64) = begin
    @inbounds for vali = x.vals
        resize!(vali, len)
    end
end

initialise!(x::Tensor; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    degree, _ = size(x)
    @inbounds for i = 1:degree
        initialise!(x.vals[i], isntt=isntt)
    end

    x.isPQ[] = isPQ
    x.auxQ[] = auxQ
end

#=================================================================================================#

"""
RLEV is a struct for RLEV ciphertext.
"""
struct RLEV
    glen::Int64
    stack::Vector{RLWE}

    function RLEV(stack::Vector{RLWE})
        @inbounds for rlwe = stack
            @assert rlwe.auxQ[] == 0 "The auxiliary modulus should be zero for RLEV encryption."
        end
        new(length(stack), stack)
    end

    function RLEV(N::Int64, len::Int64, glen::Int64; isntt::Bool=true, isPQ::Bool=false)
        new(glen, [RLWE(N, len, isntt=isntt, isPQ=isPQ) for _ = 1:glen])
    end
end

Base.:copy(x::RLEV) = RLEV(copy(x.stack))

Base.:copy!(dst::RLEV, src::RLEV) = begin
    @assert dst.glen == src.glen "The length of input and output ciphertexts should match."
    copy!(dst.stack, src.stack)
end

Base.:similar(ct::RLEV) = RLEV(ct.glen, [similar(stacki) for stacki = ct.stack])

Base.:resize!(ct::RLEV, len::Int64) = begin
    @inbounds for stacki = ct.stack
        resize!(stacki, len)
    end
end

initialise!(x::RLEV; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    @inbounds for i = 1:x.len
        initialise!(x.stack[i], isntt=isntt, isPQ=isPQ, auxQ=auxQ)
    end
end

#=================================================================================================#

"""
RGSW is a struct for RGSW ciphertext.
"""
struct RGSW
    basketb::RLEV
    basketa::RLEV

    function RGSW(basketb::RLEV, basketa::RLEV)
        @assert basketb.glen == basketa.glen "The length of input and output ciphertexts should match."
        new(basketb, basketa)
    end

    function RGSW(N::Int64, len::Int64, glen::Int64; isntt::Bool=true, isPQ::Bool=false)
        new(RLEV(N, len, glen, isntt=isntt, isPQ=isPQ), RLEV(N, len, glen, isntt=isntt, isPQ=isPQ))
    end
end

Base.:copy(x::RGSW) = RGSW(copy(x.basketb), copy(x.basketa))
Base.:copy!(dst::RGSW, src::RGSW) = begin
    copy!(dst.basketb, src.basketb)
    copy!(dst.basketa, src.basketa)
end

Base.:similar(ct::RGSW) = RGSW(similar(ct.basketb), similar(ct.basketa))

Base.:resize!(ct::RGSW, len::Int64) = begin
    resize!(ct.basketb, len)
    resize!(ct.basketa, len)
end

initialise!(x::RGSW; isntt::Bool=true, isPQ::Bool=false, auxQ::UInt64=UInt64(0)) = begin
    initialise!(x.basketb, isntt=isntt, isPQ=isPQ, auxQ=auxQ)
    initialise!(x.basketa, isntt=isntt, isPQ=isPQ, auxQ=auxQ)
end

export PlainConst, PlainPoly, PlainText, RLWE, Tensor, RLEV, RGSW, initialise!