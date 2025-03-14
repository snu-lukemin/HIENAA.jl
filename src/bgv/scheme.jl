"""
    BGVScheme(sketch::BGVParamSketch)
    BGVScheme(param::BGVParameters)
    BGVScheme(oper::BGVOperator)

BGVScheme is a struct for the BGV homomorphic encryption scheme. It contains the following fields: oper, entor, rlk, and atk.
- oper: BGVOperator
- entor: Union{Missing,Encryptor}
- rlk: Union{Missing,RLEV}
- atk: Dict{Int64,RLEV}
"""
mutable struct BGVScheme <: HEScheme
    oper::BGVOperator
    entor::Union{Missing,Encryptor}
    rlk::Union{Missing,RLEV}
    atk::Dict{Int64,RLEV}

    BGVScheme(sketch::BGVParamSketch) = BGVScheme(BGVParameters(sketch))

    function BGVScheme(param::BGVParameters)
        oper = BGVOperator(param)

        new(oper, missing, missing, Dict{Int64,RLEV}())
    end

    BGVScheme(oper::BGVOperator) = new(oper, missing, missing, Dict{Int64,RLEV}())
end

#================================================================================#
##################################### Matrix #####################################
#================================================================================#

(::Type{PlainMatrix})(M::Matrix{<:Integer}, level::Integer, oper::BGVOperator; BSGSparam::Tuple{Int64,Int64}=(0, 0))::PlainMatrix = begin
    row, col = size(M)

    if row ≠ col
        throw(DomainError("Currently, only square matrices are supported."))
    end
    if ismissing(oper.packer)
        throw(ErrorException("The operator must have a packer."))
    end
    if !isa(oper.packer, IntPackerSubring)
        throw(ErrorException("The packer must be a subring packer."))
    end
    if oper.packer.k % row ≠ 0
        throw(DomainError("The number of rows must divide the packing parameter."))
    end

    res = Dict{Int64,PlainPoly}()
    mattmp = Vector{UInt64}(undef, row)
    n1, n2 = BSGSparam == (0, 0) ? get_bsgs_param(row) : BSGSparam
    if n1 * n2 < row
        throw(DomainError("The BSGS parameters are too small."))
    end

    for i = 0:n2-1
        for j = 0:n1-1
            i * n1 + j ≥ row && break

            @inbounds for k = 0:row-1
                mattmp[k+1] = M[k+1, mod(k - i * n1 - j, row)+1]
            end
            all(iszero, mattmp) && continue

            circshift!(mattmp, -i * n1)
            res[i*n1+j] = encode(mattmp, oper, level=level, isntt_friendly=true)
        end
    end

    PlainMatrix(res, (n1, n2))
end

mul(M::PlainMatrix, x::BGV, atk::Dict{Int64,RLEV}, oper::BGVOperator; islazy::Bool=false)::BGV = begin
    res = similar(x)
    mul_to!(res, M, x, atk, oper, islazy=islazy)
    res
end

mul_to!(res::BGV, M::PlainMatrix, x::BGV, atk::Dict{Int64,RLEV}, oper::BGVOperator; islazy::Bool=false)::Nothing = begin
    # Sanity check.
    xlevel = x.level[]
    if xlevel ≤ 0
        throw(DomainError("The input ciphertext should be at least at level 1."))
    end

    now_Qlen, now_auxQ = oper.Qatlevel[xlevel+1]
    tar_Qlen, tar_auxQ = oper.Qatlevel[xlevel]
    if now_Qlen ≠ length(x) || now_auxQ ≠ x.val.auxQ[]
        throw(DomainError("The input ciphertext length does not match the parameters.."))
    end

    # Convert the input ciphertext to BFV.
    resize!(res.val, now_Qlen)
    _MSB_pack_to!(res.val, x.val, oper)

    # Scale the ciphertext to NTT-friendly modulus.
    operQ = oper.operQ
    scale_to!(res.val, res.val, now_Qlen, operQ)

    # Perform the matrix multiplication.
    _mul_RLWE!(M, res.val, atk, oper)

    if islazy
        # Scale back the ciphertext.
        scale_to!(res.val, res.val, now_Qlen, operQ, auxQ=now_auxQ)

        # Update the level.
        res.level[] = xlevel
    else
        # Rescale the ciphertext.
        scale_to!(res.val, res.val, tar_Qlen, operQ, auxQ=tar_auxQ)

        # Update the level.
        res.level[] = xlevel - 1
    end

    # Convert the output ciphertext to BGV.
    _LSB_pack_to!(res.val, res.val, oper)

    return nothing
end

#================================================================================#
############################# Polynomial Evaluation ##############################
#================================================================================#

(::Type{PolyCoeffs})(coeffs::Vector{<:Integer}, oper::BGVOperator; type::Symbol=:monomial)::PolyCoeffs = begin
    if type ≠ :monomial
        throw(DomainError("Currently BGV scheme only supports monomial basis."))
    end

    res = Dict{Int64,UInt64}()
    for i = eachindex(coeffs)
        coeffsi = Bred(coeffs[i], oper.ptxt_modulus)
        if coeffsi ≠ 0
            res[i-1] = coeffsi
        end
    end

    PolyCoeffs(res, :monomial)
end

function evaluate(poly::PolyCoeffs, x::BGV, rlk::RLEV, oper::BGVOperator)::BGV
    if poly.type ≠ :monomial
        throw(DomainError("Currently BGV scheme only supports monomial basis."))
    end

    # Compute the basis for babystep.
    basis = compute_basis(poly, x, rlk, oper)

    # Evaluate the polynomial using BSGS algorithm.
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    maxdeg = (1 << maxlevel) - 1
    res = evalrecurse(0, maxdeg, poly, basis, rlk, oper)

    if isnothing(res)
        throw(ErrorException("The polynomial evaluation failed."))
    end

    return res
end

function compute_basis(poly::PolyCoeffs, x::BGV, rlk::RLEV, oper::BGVOperator)::Dict{Int64,BGV}
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    maxdeg = (1 << maxlevel) - 1
    babydeg = 1 << ceil(Int64, maxlevel / 2)
    babylevel = trailing_zeros(babydeg)
    basis = Dict{Int64,BGV}()

    # Compute the power-of-two monomials.
    basis[1] = deepcopy(x)
    for i = 1:maxlevel-1
        basis[1<<i] = mul(basis[1<<(i-1)], basis[1<<(i-1)], rlk, oper)
    end

    @inbounds for i = 1:babylevel-1
        # Compute the power-of-two monomial.
        for j = 1:(1<<i)-1
            # Checks if the polynomial is sparse.
            checkj = !isassigned(poly, (1 << i) + j)

            # Used in baby step generation?
            for idx = i+1:babylevel-1
                checkj = checkj && !isassigned(poly, (1 << idx) + (1 << i) + j)
            end
            # Used in the final evaluation?
            for idx = 0:maxdeg÷babydeg-2
                checkj = checkj && !isassigned(poly, idx * babydeg + j)
            end
            checkj && continue

            # Compute the monomial if needed in the future.
            basis[(1<<i)+j] = mul(basis[j], basis[1<<i], rlk, oper)
        end
    end

    return basis
end

function evalrecurse(lo::Int64, hi::Int64, poly::PolyCoeffs, basis::Dict{Int64,BGV}, rlk::RLEV, oper::BGVOperator)::Union{BGV,Nothing}
    # Hyperparameters.
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    maxdeg = (1 << maxlevel) - 1
    babydeg = 1 << ceil(Int64, maxlevel / 2)
    babylevel = trailing_zeros(babydeg)

    # Current degree.
    d = hi - lo + 1
    if !ispow2(d)
        throw(DomainError("The degree should be a power of 2."))
    end

    # Corner case. (Does anybody really want to evaluate a constant polynomial?)
    if d == 1
        res = similar(basis[1])
        isassigned(poly, lo) && add_to!(res, res, poly[lo], oper)
        return res
    end

    # Temporary variables.
    ctlen = length(basis[1])
    tmp = oper.ct_buff[3][1:ctlen]

    if d == babydeg
        lo > deg && return nothing
        coeffi = PlainConst(ctlen)

        res = similar(basis[1])
        if hi == maxdeg
            # Babystep computation, at the highest degree.
            # We need to perform BSGS to the smallest degree to minimise the level consumption.
            for i = 0:babylevel
                d, halfd = 1 << i, 1 << (i - 1)
                maxdeg - d ≥ deg && continue

                if i == 0
                    encode_to!(coeffi, poly[maxdeg-d+1], oper, level=res.level[])
                    mul_to!(res, coeffi, basis[1], oper)
                else
                    if isassigned(poly, maxdeg - d + 1)
                        encode_to!(coeffi, poly[maxdeg-d+1], oper, level=res.level[])
                        add_to!(res, res, coeffi, oper)
                    end

                    for j = 1:min(halfd - 1, deg + d - maxdeg - 1)
                        if !isassigned(poly, maxdeg - d + j + 1)
                            continue
                        else
                            encode_to!(coeffi, poly[maxdeg-d+j+1], oper, level=basis[j].level[])
                            mul_to!(tmp, coeffi, basis[j], oper)
                            add_to!(res, res, tmp, oper)
                        end
                    end

                    i == babylevel && break
                    mul_to!(res, res, basis[1<<i], rlk, oper)
                end
            end
        else
            # Babystep computation.
            if isassigned(poly, lo)
                encode_to!(coeffi, poly[lo], oper, level=res.level[])
                add_to!(res, res, coeffi, oper)
            end
            for i = 1:min(babydeg - 1, deg - lo)
                if isassigned(poly, lo + i)
                    encode_to!(coeffi, poly[lo+i], oper, level=basis[i].level[])
                    mul_to!(tmp, coeffi, basis[i], oper)
                    add_to!(res, res, tmp, oper)
                end
            end
        end

        return res
    else
        # Giantstep computation.
        halfd = d >> 1

        lo_res = evalrecurse(lo, lo + halfd - 1, poly, basis, rlk, oper)
        hi_res = evalrecurse(lo + halfd, hi, poly, basis, rlk, oper)

        if isnothing(hi_res)
            isnothing(lo_res) && return nothing
        else
            mul_to!(tmp, hi_res, basis[halfd], rlk, oper)
            add_to!(lo_res, lo_res, tmp, oper)
        end

        return lo_res
    end
end