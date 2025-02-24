
get_bsgs_param(n::Int64) = begin
    n1 = ceil(Int64, sqrt(n))
    n2 = ceil(Int64, n / n1)
    n1, n2
end

"""
    PlainMatrix(M::Matrix{UInt64}, scheme::IntScheme; BSGSparam::Tuple{Int64,Int64}=(0, 0))
    PlainMatrix(M::Matrix{<:Complex{<:AbstractFloat}}, scheme::CKKSScheme; BSGSparam::Tuple{Int64,Int64}=(0, 0))

A struct for diagonal-packed matrices.
The input matrix should be a square matrix, and the number of rows should divide the packing parameter.
"""
struct PlainMatrix
    diag::Dict{Int64,PlainPoly}
    BSGSparam::Tuple{Int64,Int64}

    function PlainMatrix(M::Matrix{UInt64}, scheme::IntScheme; BSGSparam::Tuple{Int64,Int64}=(0, 0))
        oper = scheme.oper
        row, col = size(M)
        @assert row == col "Currently, only square matrices are supported."
        @assert !ismissing(oper.packer) "The operator must have a packer."
        @assert typeof(oper.packer) == IntPackerSubring "The packer must be a subring packer."
        @assert oper.packer.k % row == 0 "The number of rows must divide the packing parameter."

        res = Dict{Int64,PlainPoly}()
        mattmp = Vector{UInt64}(undef, row)
        n1, n2 = BSGSparam == (0, 0) ? get_bsgs_param(row) : BSGSparam
        @assert n1 * n2 ≥ row "Something is wrong with the BSGS parameters."

        for i = 0:n2-1
            for j = 0:n1-1
                i * n1 + j ≥ row && break

                @inbounds for k = 0:row-1
                    mattmp[k+1] = M[k+1, mod(k - i * n1 - j, row)+1]
                end
                all(iszero, mattmp) && continue

                circshift!(mattmp, -i * n1)
                res[i*n1+j] = encode(mattmp, scheme)
            end
        end

        new(res, (n1, n2))
    end

    function PlainMatrix(M::Matrix{<:Complex{<:AbstractFloat}}, scheme::CKKSScheme; BSGSparam::Tuple{Int64,Int64}=(0, 0))
        oper = scheme.oper
        row, col = size(M)
        @assert row == col "Currently, only square matrices are supported."
        @assert !ismissing(oper.packer) "The operator must have a packer."
        @assert typeof(oper.packer) == ComplexPackerPow2 "The packer must be a subring packer."
        @assert oper.packer.k % row == 0 "The number of rows must divide the packing parameter."

        res = Dict{Int64,PlainPoly}()
        mattmp = Vector{ComplexDF64}(undef, row)
        n1, n2 = BSGSparam == (0, 0) ? get_bsgs_param(row) : BSGSparam
        @assert n1 * n2 ≥ row "Something is wrong with the BSGS parameters."

        for i = 0:n2-1
            for j = 0:n1-1
                i * n1 + j ≥ row && break

                @inbounds for k = 0:row-1
                    mattmp[k+1] = M[k+1, mod(k - i * n1 - j, row)+1]
                end
                all(x -> isapprox(x, 0, atol=1 / oper.scaling_factor), mattmp) && continue

                circshift!(mattmp, -i * n1)
                res[i*n1+j] = encode(mattmp, scheme)
            end
        end

        new(res, (n1, n2))
    end
end

"""
    get_required_key_list(M::PlainMatrix)

Get the list of keys required for the multiplication with the given matrix.
"""
function get_required_key_list(M::PlainMatrix)
    res = Vector{Tuple{Int64}}(undef, 0)
    n1, n2 = M.BSGSparam
    for i = 0:n2-1
        check = false
        for j = 1:n1-1
            if haskey(M.diag, i * n1 + j)
                check = true
                if (j,) ∉ res
                    push!(res, (j,))
                end
            end
        end
        if check && i ≠ 0
            push!(res, (i * n1,))
        end
    end

    sort!(res, by=x -> x[1])

    res
end

"""
    mul(M::PlainMatrix, x::HECiphertext, scheme::HEScheme)

Multiply the given matrix with the ciphertext.
"""
function mul(M::PlainMatrix, x::HECiphertext, scheme::HEScheme)
    oper = scheme.oper

    res = similar(x)
    tmp = oper.ct_buff[3][1:length(x.val)]
    tmp2 = oper.ct_buff[4][1:length(x.val)]
    tmp.level[], tmp2.level[] = x.level[], x.level[]

    adec = decompose_a(x, oper)
    n1, n2 = M.BSGSparam

    rotx = Vector{typeof(x)}(undef, n1)
    rotx[1] = x

    # TODO double hoisting.

    diag = M.diag
    for i = 0:n2-1
        initialise!(tmp2.val)
        check = false
        for j = 0:n1-1
            !haskey(diag, i * n1 + j) && continue
            check = true

            if !isassigned(rotx, j + 1)
                rotx[j+1] = hoisted_rotate(adec, x, (j,), scheme)
            end

            if j == 0
                mul_to!(tmp2, diag[i*n1+j], rotx[j+1], oper, islazy=true)
            else
                mul_to!(tmp, diag[i*n1+j], rotx[j+1], oper, islazy=true)
                add_to!(tmp2, tmp, tmp2, oper)
            end
        end

        if i == 0
            copy!(res, tmp2)
        else
            !check && continue
            rotate_to!(tmp, tmp2, (i * n1,), scheme)
            add_to!(res, tmp, res, oper)
        end
    end

    rescale_to!(res, res, oper)

    res
end

export PlainMatrix, get_required_key_list