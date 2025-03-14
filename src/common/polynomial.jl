"""
    PolyCoeffs(coeffs::Vector{<:Number}, scheme::HEScheme; type::Symbol=:monomial)

A polynomial representation in the form of a dictionary of coefficients. 
The polynomial can be of two types: monomial or Chebyshev. The monomial basis is the standard basis for polynomials, while the Chebyshev basis is a basis of orthogonal polynomials. 
The Chebyshev basis is often used in the CKKS scheme for polynomial evaluation with large degree.
"""
struct PolyCoeffs
    coeffs::Dict{Int64, <:Number}
    type::Symbol # :monomial, :chebyshev.

    PolyCoeffs(coeffs::Dict{Int64, <:Number}, type::Symbol=:monomial)::PolyCoeffs = new(coeffs, type)
    
    PolyCoeffs(coeffs::Vector{<:Number}, scheme::HEScheme; type::Symbol=:monomial)::PolyCoeffs = 
        PolyCoeffs(coeffs, scheme.oper, type=type)
end

"""
    degree(coeffs::PolyCoeffs)

Returns the degree of the polynomial.
"""
degree(coeffs::PolyCoeffs)::Int64 = maximum(keys(coeffs.coeffs))
Base.:isassigned(coeffs::PolyCoeffs, i::Integer)::Bool = haskey(coeffs.coeffs, i)
Base.:getindex(coeffs::PolyCoeffs, i::Integer)::Number = getindex(coeffs.coeffs, i)

"""
    evaluate(coeffs::PolyCoeffs, x::HECiphertext, scheme::HEScheme)
    evaluate(coeffs::Vector{PolyCoeffs}, x::HECiphertext, scheme::HEScheme)

Evaluates the polynomial(s) `coeffs` at the ciphertext `x` in depth-optimal manner.
"""
function evaluate(poly::PolyCoeffs, x::HECiphertext, scheme::HEScheme)::HECiphertext
    rlk, oper = scheme.rlk, scheme.oper
    if ismissing(rlk)
        throw(ErrorException("The relinearisation key is required for polynomial evaluation."))
    end
    evaluate(poly, x, rlk, oper)
end