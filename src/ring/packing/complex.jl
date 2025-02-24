struct ComplexPackerPow2
    m::Int64
    N::Int64
    d::Int64
    k::Int64
    cube::Vector{Int64}
    cubegen::Vector{Int64}
    ζ::Vector{ComplexDF64}
    group::Vector{Int64}
    buff::Vector{ComplexDF64}

    function ComplexPackerPow2(m::Int64)
        N = m >> 1
        d, k = 1, N >> 1
        cube = [N >> 1]
        cubegen = [5]

        setprecision(64 * 3)
        ζ = ComplexDF64.(exp.(2big(pi) * im .* (0:2N-1) / 2N))
        group = [powermod(5, (N >> 1) - i, 2N) for i = 0:N-1]
        buff = zeros(ComplexDF64, N >> 1)
        new(m, N, d, k, cube, cubegen, ζ, group, buff)
    end
end

function pack_to!(res::AbstractVector{Int128}, msg::AbstractVector{<:Complex{<:AbstractFloat}}, Δ::Real, packer::ComplexPackerPow2)
    N, ζ, group, buff = packer.N, packer.ζ, packer.group, packer.buff
    n = length(msg)
    @assert N % 2n == 0

    @views @. buff[1:n] = msg * Δ / n
    @inbounds for idx = reverse(0:trailing_zeros(n))
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
        res[i*N÷2n+1] = Int128(round(real(buff[i+1])))
        res[(i+n)*N÷2n+1] = Int128(round(imag(buff[i+1])))
    end
end

function unpack_to!(res::AbstractVector{ComplexDF64}, pt::AbstractVector{Int128}, Δ::Integer, packer::ComplexPackerPow2)
    N, ζ, group, buff = packer.N, packer.ζ, packer.group, packer.buff

    @assert N == length(pt)
    n = N >> 1

    @inbounds for i = 1:n
        buff[i] = ComplexDF64(pt[i], pt[i+n]) / Δ
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
    @views @. res = buff[1:reslen]
end

const ComplexPacker = Union{ComplexPackerPow2}

(::Type{ComplexPacker})(param::RingParam) = begin
    @assert typeof(param) ≠ CyclicParam "Cyclic parameters are not supported."

    m = param.m

    if ispow2(m)
        ComplexPackerPow2(m)
    elseif typeof(param) == SubringParam
        @error "Unsupported ring parameters."
    else
        @error "Unsupported ring parameters."
    end
end