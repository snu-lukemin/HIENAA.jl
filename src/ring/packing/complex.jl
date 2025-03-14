abstract type ComplexPacker end

struct ComplexPackerPow2 <: ComplexPacker
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    ζ::Vector{ComplexBF}
    group::Vector{Int64}
    buff::Vector{ComplexBF}

    function ComplexPackerPow2(m::Int64)::ComplexPackerPow2
        N = m >> 1
        d, k = 1, N >> 1
        cube = [N >> 1]
        cubegen = [5]

        setprecision(192)
        ζ = exp.(2big(pi) * im .* (0:2N-1) / 2N)
        group = [powermod(5, (N >> 1) - i, 2N) for i = 0:N-1]
        buff = zeros(ComplexBF, N >> 1)
        new(m, N, d, k, cube, cubegen, ζ, group, buff)
    end
end

function pack_to!(res::AbstractVector{Int128}, msg::AbstractVector{<:Union{Real,Complex{<:Real}}}, Δ::Real, packer::ComplexPackerPow2)::Nothing
    N, ζ, group, buff = packer.N, packer.ζ, packer.group, packer.buff
    n = length(msg)
    if N % 2n ≠ 0
        throw(DimensionMismatch("The length of the message should be a divisor of the packing parameter."))
    end

    @views @. buff[1:n] = msg * Δ / n
    @inbounds for idx = reverse(1:trailing_zeros(n))
        len = 1 << idx
        for i = 0:len:n-1
            lenh, lenQ = len >> 1, len << 2
            gap = 2N ÷ lenQ
            for j in 0:lenh-1
                idx1, idx2 = i + j + 1, i + j + lenh + 1
                ζi = ζ[(lenQ-(group[j+1]&(lenQ-1)))*gap+1]
                buff[idx1], buff[idx2] = buff[idx1] + buff[idx2], (buff[idx1] - buff[idx2]) * ζi
            end
        end
    end
    @views scramble!(buff[1:n], 2)

    @. res = 0
    @inbounds for i = 0:n-1
        res[i*N÷2n+1] = round(Int128, real(buff[i+1]))
        res[(i+n)*N÷2n+1] = round(Int128, imag(buff[i+1]))
    end

    return nothing
end

function unpack_to!(res::AbstractVector{ComplexBF}, pt::AbstractVector{Int128}, Δ::Real, packer::ComplexPackerPow2)::Nothing
    N, ζ, group, buff = packer.N, packer.ζ, packer.group, packer.buff

    if N ≠ length(pt)
        throw(DimensionMismatch("The length of the packed vector should be the same as the packing parameter."))
    end
    n = N >> 1

    for i = 1:n
        buff[i] = pt[i] + im * pt[i+n]
    end

    @views scramble!(buff[1:n], 2)
    @inbounds for idx = 1:trailing_zeros(n)
        len = 1 << idx
        for i in 0:len:n-1
            lenh, lenQ = len >> 1, len << 2
            gap = 2N ÷ lenQ
            for j in 0:lenh-1
                tmp = ζ[(group[j+1]&(lenQ-1))*gap+1] * buff[i+j+lenh+1]
                buff[i+j+1], buff[i+j+lenh+1] = buff[i+j+1] + tmp, buff[i+j+1] - tmp
            end
        end
    end

    reslen = length(res)
    @views @. res = buff[1:reslen] / Δ

    return nothing
end

struct ComplexPackerArb <: ComplexPacker
    function ComplexPackerArb(m::Int64)::ComplexPackerArb
        error("Not implemented yet.")
    end
end

struct ComplexPackerSubring <: ComplexPacker
    function ComplexPackerSubring(m::Int64, d::Int64)::ComplexPackerSubring
        error("Not implemented yet.")
    end
end

(::Type{ComplexPacker})(param::RingParam)::ComplexPacker = begin
    if isa(param, CyclicParam)
        throw(DomainError("Cyclic parameters are not supported."))
    end

    m = param.m

    if ispow2(m)
        ComplexPackerPow2(m)
    elseif isa(param, SubringParam)
        ComplexPackerSubring(m, param.d)
    else
        ComplexPackerArb(m)
    end
end