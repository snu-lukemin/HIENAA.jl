get_bsgs_param(n::Int64)::Tuple{Int64,Int64} = begin
    n1 = ceil(Int64, sqrt(n))
    n2 = ceil(Int64, n / n1)
    n1, n2
end

"""
    PlainMatrix(M::Matrix{<:Number}, level::Integer, scheme::HEScheme; BSGSparam::Tuple{Int64,Int64}=(0, 0))

A struct for diagonal-packed matrices.
The input matrix should be a square matrix, and the number of rows should divide the packing parameter.
We require the evaluated level for the best performance / noise.
"""
struct PlainMatrix
    diag::Dict{Int64,PlainPoly}
    BSGSparam::Tuple{Int64,Int64}

    PlainMatrix(diag::Dict{Int64,PlainPoly}, BSGSparam::Tuple{Int64,Int64})::PlainMatrix = new(diag, BSGSparam)

    PlainMatrix(M::Matrix{<:Number}, level::Integer, scheme::HEScheme; BSGSparam::Tuple{Int64,Int64}=(0, 0))::PlainMatrix = 
        PlainMatrix(M, level, scheme.oper, BSGSparam=BSGSparam)
end

"""
    get_required_key_list(M::PlainMatrix)

Get the list of keys required for the multiplication with the given matrix.
"""
function get_required_key_list(M::PlainMatrix)::Vector{Tuple{Int64}}
    res = Vector{Tuple{Int64}}(undef, 0)
    n1, n2 = M.BSGSparam
    for i = 0:n2-1
        check = false
        for j = 0:n1-1
            if haskey(M.diag, i * n1 + j)
                check = true
                if j ≠ 0 && (j,) ∉ res
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
function mul(M::PlainMatrix, x::HECiphertext, scheme::HEScheme)::HECiphertext
    atk, oper = scheme.atk, scheme.oper
    mul(M, x, atk, oper)
end

function _mul_RLWE!(M::PlainMatrix, input::RLWE, atk::Dict{Int64,RLEV}, oper::HEOperator)::Nothing
    if input.auxQ[] ≠ 0
        throw(DomainError("The input ciphertext should be NTT-friendly."))
    end

    # Define structs.
    operQ = oper.operQ
    packer = oper.packer

    # Decompose a part.
    adec = decompose(PlainPoly(input.a), operQ)
    n1, n2 = M.BSGSparam

    # Generate the rotation vectors.
    rotx = Vector{RLWE}(undef, n1)
    rotx[1] = deepcopy(input)

    # Generate the output ciphertexts.
    ctlen = length(input)
    tmp = oper.ct_buff[3].val[1:ctlen]
    tmp2 = oper.ct_buff[4].val[1:ctlen]

    # Perform BSGS algorithm.
    cube, cubegen, m = packer.cube, packer.cubegen, packer.m
    diag = M.diag
    for i = 0:n2-1
        initialise!(tmp2)
        checkj = false
        for j = 0:n1-1
            !haskey(diag, i * n1 + j) && continue

            autidx = gen_power_modmul(cubegen, cube, (j,), m)
            !isassigned(rotx, j + 1) && (rotx[j+1] = hoisted_automorphism(adec, input, autidx, atk[autidx], operQ))

            if !checkj
                mul_to!(tmp2, diag[i*n1+j], rotx[j+1], operQ)
                checkj = true
            else
                mul_to!(tmp, diag[i*n1+j], rotx[j+1], operQ)
                add_to!(tmp2, tmp, tmp2, operQ)
            end
        end

        !checkj && continue

        if i == 0
            copy!(input, tmp2)
        else
            autidx = gen_power_modmul(cubegen, cube, (i * n1,), m)
            automorphism_to!(tmp, tmp2, autidx, atk[autidx], operQ)
            add_to!(input, tmp, input, operQ)
        end
    end

    return nothing
end