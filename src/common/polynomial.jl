"""
    PolyCoeffs(coeffs::Vector{UInt64}, scheme::IntScheme; type::Symbol=:monomial)
    PolyCoeffs(coeffs::Vector{<:AbstractFloat}, scheme::CKKSScheme; type::Symbol=:monomial)

A polynomial representation in the form of a vector of coefficients. 
The coefficients are encoded as `PlainConst` elements, which are the plaintext representation of the coefficients. 
The polynomial can be of two types: monomial or Chebyshev. The monomial basis is the standard basis for polynomials, while the Chebyshev basis is a basis of orthogonal polynomials. 
The Chebyshev basis is often used in the CKKS scheme for polynomial evaluation with large degree.
"""
struct PolyCoeffs
    coeffs::Vector{PlainConst} # CKKS offers evaluation of complex numbers, but is not really utilised.
    type::Symbol # :monomial, :chebyshev.

    # For Int HE schemes.
    function PolyCoeffs(coeffs::Vector{UInt64}, scheme::IntScheme; type::Symbol=:monomial)
        @assert type == :monomial "Currently Int HE schemes only support monomial basis."
        res = Vector{PlainConst}(undef, length(coeffs))
        for i = eachindex(res)
            if _Bred(coeffs[i], scheme.oper.ptxt_modulus) ≠ 0
                res[i] = encode(coeffs[i], scheme)
            end
        end

        cnt = length(res)
        while true
            if !isassigned(res, cnt)
                cnt -= 1
            else
                break
            end
        end
        resize!(res, cnt)

        new(res, :monomial)
    end

    # For CKKS polynomials.
    function PolyCoeffs(coeffs::Vector{<:AbstractFloat}, scheme::CKKSScheme; type::Symbol=:monomial)
        @assert type == :monomial || type == :chebyshev "Unknown polynomial types."
        res = Vector{PlainConst}(undef, length(coeffs))
        for i = eachindex(res)
            if abs(coeffs[i]) ≥ 1 / scheme.oper.scaling_factor
                res[i] = encode(coeffs[i], scheme)
            end
        end

        cnt = length(res)
        while true
            if !isassigned(res, cnt)
                cnt -= 1
            else
                break
            end
        end
        resize!(res, cnt)

        new(res, type)
    end
end

"""
    degree(coeffs::PolyCoeffs)

Returns the degree of the polynomial.
"""
degree(coeffs::PolyCoeffs) = length(coeffs.coeffs) - 1
Base.:isassigned(coeffs::PolyCoeffs, i::Integer) = isassigned(coeffs.coeffs, i + 1)
Base.:getindex(coeffs::PolyCoeffs, i::Integer) = getindex(coeffs.coeffs, i + 1)

"""
    evaluate(coeffs::PolyCoeffs, x::HECiphertext, scheme::HEScheme)
    evaluate(coeffs::Vector{PolyCoeffs}, x::HECiphertext, scheme::HEScheme)

Evaluates the polynomial(s) `coeffs` at the ciphertext `x` in depth-optimal manner.
"""
function evaluate(poly::PolyCoeffs, x::HECiphertext, scheme::HEScheme)
    # Compute the maximum level.
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    halfdeg = 1 << (maxlevel - 1)
    first_half = Vector{typeof(x)}(undef, halfdeg)

    # Compute the monomials for the first half.
    first_half[1] = x
    @inbounds for i = 1:maxlevel-2
        # Compute the power-of-two monomial.
        first_half[1<<i] = mul(first_half[1<<(i-1)], first_half[1<<(i-1)], scheme)
        for j = 1:(1<<i)-1
            # Checks if the polynomial is sparse.
            checkj = !isassigned(poly, (1 << i) + j)
            for idx = i+1:maxlevel-2
                checkj = checkj && !isassigned(poly, (1 << idx) + (1 << i) + j)
            end
            if j + halfdeg ≤ deg
                checkj = checkj && !isassigned(poly, j + halfdeg)
            end
            checkj && continue

            # Compute the monomial if needed in the future.
            first_half[(1<<i)+j] = mul(first_half[j], first_half[1<<i], scheme)
        end
    end
    first_half[1<<(maxlevel-1)] = mul(first_half[1<<(maxlevel-2)], first_half[1<<(maxlevel-2)], scheme)

    _poly_eval_from_monomial(poly, first_half, x, scheme)
end

function evaluate(poly::Vector{PolyCoeffs}, x::HECiphertext, scheme::HEScheme)
    maxdeg = maximum([degree(polyi) for polyi = poly]) - 1
    maxlevel = trailing_zeros(Base._nextpow2(maxdeg + 1))
    halfdeg = 1 << (maxlevel - 1)
    first_half = Vector{typeof(x)}(undef, halfdeg)

    # Compute the monomials for the first half.
    first_half[1] = x
    @inbounds for i = 1:maxlevel-2
        # Compute the power-of-two monomial.
        first_half[1<<i] = mul(first_half[1<<(i-1)], first_half[1<<(i-1)], scheme)
        for j = 1:(1<<i)-1
            # Checks if the polynomial is sparse.
            checkj = true
            for idx1 = eachindex(poly)
                (1 << i) + j > degree(poly[idx1]) && continue
                checkj = checkj && !isassigned(poly[idx1], (1 << i) + j)
                for idx2 = i+1:maxlevel-2
                    (1 << idx2) + (1 << i) + j > degree(poly[idx1]) && continue
                    checkj = checkj && !isassigned(poly[idx1], (1<<idx2)+(1<<i)+j)
                end
                if j + halfdeg ≤ deg
                    checkj = checkj && !isassigned(poly[idx1], j+halfdeg)
                end
            end
            checkj && continue

            # Compute the monomial if needed in the future.
            first_half[(1<<i)+j] = mul(first_half[j], first_half[1<<i], scheme)
        end
    end
    first_half[1<<(maxlevel-1)] = mul(first_half[1<<(maxlevel-2)], first_half[1<<(maxlevel-2)], scheme)

    res = Vector{typeof(x)}(undef, length(poly))
    for i = eachindex(poly)
        deg = degree(poly[i])
        maxleveli = trailing_zeros(Base._nextpow2(deg + 1))
        halfdegi = 1 << (maxleveli - 1)

        @views res[i] = _poly_eval_from_monomial(poly[i], first_half[1:halfdegi], x, scheme)
    end

    res
end

# Depth-optimal naïve polynomial evaluation.
function _poly_eval_from_monomial(poly::PolyCoeffs, first_half::AbstractVector{T}, x::T, scheme::HEScheme) where {T<:HECiphertext}
    deg = degree(poly)
    maxlevel = trailing_zeros(Base._nextpow2(deg + 1))
    halfdeg = 1 << (maxlevel - 1)
    @assert length(first_half) == halfdeg "first_half must have length 2^(maxlevel-1)."

    res = similar(x)
    tmp = similar(x)

    # Compute the monomials for the last half and sum, using the binary tree.
    halfhalfdeg = halfdeg >> 1
    @inbounds for i = 1:deg-halfdeg
        !isassigned(poly, i + halfdeg) && continue

        tmpi = tmp[1:length(tmp.val)]

        if i ≤ halfhalfdeg
            mul_to!(tmpi, poly[i+halfdeg], first_half[i], scheme)
        else
            checki = false
            for j = 0:maxlevel-2
                if (i >> j) & 1 == 1
                    if !checki
                        mul_to!(tmpi, first_half[1<<j], poly[i+halfdeg], scheme)
                        checki = true
                    else
                        mul_to!(tmpi, tmpi, first_half[1<<j], scheme)
                    end
                end
            end
        end

        add_to!(res, tmpi, res, scheme)
    end

    # Multiply x^halfdeg for a better computation complexity.
    mul_to!(res, res, first_half[halfdeg], scheme)

    # Compute the linear sum.
    @inbounds for i = 1:halfdeg
        !isassigned(poly, i) && continue

        mul_to!(tmp, poly[i], first_half[i], scheme)
        add_to!(res, tmp, res, scheme)
    end

    if isassigned(poly, 0)
        add_to!(res, poly[0], res, scheme)
    end

    res
end

export PolyCoeffs, degree, evaluate, evaluate_many