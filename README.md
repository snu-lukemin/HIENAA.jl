# HIENAA.jl

HIENAA (HE Implementation for Encrypted Numbers Arithmetic and Algorithms) is an open-source Julia implementation for Fully Homomorphic Encryption (FHE) schemes.

<p align="center">
  <img src="SPOTTY.png" width="350"/>
  <br> Spotty, mascot of HIENAA.jl.
</p>

**HIENAA** is a pure Julia implementation of the [[BGV](https://eprint.iacr.org/2011/277)], [[BFV](https://eprint.iacr.org/2012/144)] and [[CKKS](https://eprint.iacr.org/2016/421)] scheme. Currently, it provides 
- Ptxt-Ctxt/Ctxt-Ctxt addition/subtraction/multiplication
- (Almost) arbitrary cyclotomic ring of integers as a base ring (only for BGV/BFV)
  - The cyclotomic degree should be a multiple of at most four integers for optimisations.
- Integer packing in SIMD manner
- Subring Homomorphic Encryption [[AH17](https://eprint.iacr.org/2017/066)]
- Matrix-vector multiplication [[HS14](https://eprint.iacr.org/2014/106)] and polynomial evaluation.

## Installation
HIENAA can be installed by typing
```
] add https://github.com/snu-lukemin/HIENAA.jl
```
into the Julia REPL.

## Examples
### Encryption
```julia
# 128-bit security parameter set.
m, hw, logP, logQ = 4369, 2048, 61, 40

ring_param = CyclotomicParam(m) # Generate ring parameters.
sketch = BGVParamSketch(ring_param, logP, logQ, 2^5) # Use sketch for an easier implementation.
scheme = BGVScheme(sketch) # Compile the sketch into parameters.

# Generate secret key.
us = UniformSampler() # Define uniform sampler.
key = ternary_ringkey(us, ring_param.N, hw) # Sample sparse ternary with Hamming weight hw.
set_encryptor!(key, scheme) # Set the scheme parameters.

# Encrypt 
t, packlen = UInt64(sketch.ptxt_modulus), scheme.oper.packer.k
msg = rand(0:t-1, packlen)
pt = encode(msg, scheme) # Encode the message into a plaintext.
ct = encrypt(pt, scheme) # Encrypt the plaintext.

# Decrypt
decrypt_to!(pt, ct, scheme) # Decrypt the ciphertext
out = decode(pt, scheme) # Decode the plaintext.

# Output
println(msg[1:10] .% Int64)
println(out[1:10] .% Int64)
```

### Polynomial Evaluation
```julia
# 128-bit security parameter set.
m, hw, logP, logQ = 21845, 128, 62, 275

ring_param = CyclotomicParam(m) # Generate ring parameters.
sketch = BFVParamSketch(ring_param, logP, logQ, 2^8) # Use sketch for an easier implementation.
scheme = BFVScheme(sketch) # Compile the sketch into parameters.

# Generate secret key.
us = UniformSampler() # Define uniform sampler.
key = ternary_ringkey(us, ring_param.N, hw) # Sample sparse ternary with Hamming weight hw.
set_encryptor!(key, scheme) # Set the scheme parameters.

# Generate and set relin key.
rlk = relin_keygen(scheme)
set_relinkey!(rlk, scheme)

# Generate a polynomial.
# This is a digit extraction polynomial modulo 2^8.
coeffs = PolyCoeffs(UInt64[0, 0, 0, 0, 0, 0, 244, 0, 13], scheme)

# Encrypt
t, packlen = UInt64(sketch.ptxt_modulus), scheme.oper.packer.k
msg = rand(0:t-1, packlen)
pt = encode(msg, scheme) # Encode the message into a plaintext.
ct = encrypt(pt, scheme) # Encrypt the plaintext.

# Evaluate Polynomial.
res = evaluate(coeffs, ct, scheme)

# Decrypt
decrypt_to!(pt, res, scheme) # Decrypt the ciphertext
out = decode(pt, scheme) # Decode the plaintext.

# Output
println(msg[1:10] .% Int64)
println(out[1:10] .% Int64) # out = msg (mod 2).
```

### Matrix Multiplication
```julia
# 128-bit security parameter set.
m, hw, logP, logQ = 1 << 14, 256, 62, 150

ring_param = CyclotomicParam(m) # Generate ring parameters.
sketch = CKKSParamSketch(ring_param, logP, logQ, 1 << 40) # Use sketch for an easier implementation.
scheme = CKKSScheme(sketch) # Compile the sketch into parameters.

# Generate secret key.
us = UniformSampler() # Define uniform sampler.
key = ternary_ringkey(us, ring_param.N, hw) # Sample sparse ternary with Hamming weight hw.
set_encryptor!(key, scheme) # Set the scheme parameters.

# Prepare the matrix.
packlen = scheme.oper.packer.k # Number of the slots.
matrix = rand(ComplexF64, packlen, packlen)
M = PlainMatrix(matrix, scheme) # Encode the matrix.

# Generate rotation keys.
list = get_required_key_list(M) # Get the list of keys required for the matrix multiplication.
rtk = rotate_keygen(list, scheme) # Generate rotation keys from the list.
set_rotate_key!(list, rtk, scheme) # Set the rotation keys to the scheme.

# Encrypt
msg = rand(ComplexF64, packlen)
pt = encode(msg, scheme) # Encode the message into a plaintext.
ct = encrypt(pt, scheme) # Encrypt the plaintext.

# Matrix multiplication.
res = mul(M, ct, scheme)

# Decrypt
decrypt_to!(pt, res, scheme) # Decrypt the ciphertext
out = decode(pt, scheme) # Decode the plaintext.

# Output
println((matrix * msg)[1:10])
println(ComplexF64.(out[1:10]))
```

## License
HIENAA is licensed under the Apache 2.0 License.

## Citing
To cite HIENAA, please use the following BibTeX entry:
```tex
@misc{HIENAA.jl,
  title={{HIENAA.jl}},
  author={Seonhong Min},
  year={2025},
  howpublished = {Online: \url{https://github.com/snu-lukemin/HIENAA.jl}},
}
```

## About Spotty
Spotty is the only cub remaining from the HIENAA species, which is a special hyena with three dots on their rump. You can help the species of HIENAA to be preserved and spread by using the HIENAA.jl library. Spotty was designed by [Hyeonji Park](https://www.instagram.com/mlgng2010/). 

## References
- Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography (https://eprint.iacr.org/2016/504)
- Subring Homomorphic Encryption (https://eprint.iacr.org/2017/066)
- Efficient reductions in cyclotomic rings (https://eprint.iacr.org/2017/748)
- An Improved RNS Variant of the BFV Homomorphic Encryption Scheme (https://eprint.iacr.org/2018/117)
- Efficient Bootstrapping for Approximate Homomorphic Encryption with Non-Sparse Keys (https://eprint.iacr.org/2020/1203)
- On Polynomial Functions Modulo p^e and Faster Bootstrapping for Homomorphic Encryption (https://eprint.iacr.org/2022/1364)