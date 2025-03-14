"""
    HEBootKey

Abstract type for FHE bootstrapping key.
"""
abstract type HEBootKey end

"""
    HEBootstrapper

Abstract type for FHE bootstrapping.
"""
abstract type HEBootstrapper end

"""
    HEBootParamSketch

Abstract type for FHE bootstrapping parameter sketch.
"""
abstract type HEBootParamSketch <: HEParamSketch end

"""
    HEBootParameters

Abstract type for FHE bootstrapping parameters.
"""
abstract type HEBootParameters <: HEParameters end

"""
    bootstrap_keygen(bootparam::HEBootParameters, entor::SKEncryptor)

Generates a bootstrapping key for the homomorphic encryption scheme `scheme` with the bootstrapping parameters `bootparam` and the secret key encryptor `entor`.
"""
bootstrap_keygen(bootparam::HEBootParameters, scheme::HEScheme)::HEBootKey = begin
    oper, entor = scheme.oper, scheme.entor
    if ismissing(entor)
        throw(ErrorException("The encryptor is not generated."))
    end
    if !isa(entor, SKEncryptor)
        throw(ErrorException("The Encryptor should be a secret key encryptor."))
    end

    bootstrap_keygen(bootparam, oper, entor)
end

"""
   set_bootstrapper!(bootkey::HEBootKey, scheme::HEScheme)

Sets the bootstrapper for the scheme.
"""
set_bootstrapper!(bootkey::HEBootKey, scheme::HEScheme)::Nothing = begin
    set_automorphism_key!(collect(keys(bootkey.atk)), collect(values(bootkey.atk)), scheme)
    set_relinkey!(bootkey.rlk, scheme)
    scheme.boot = HEBootstrapper(bootkey, scheme.oper)

    return nothing
end

"""
    bootstrap(ct::HECiphertext, scheme::HEScheme)

Bootstraps a ciphertext `ct` using the bootstrapper `boot`.
"""
bootstrap(ct::HECiphertext, scheme::HEScheme)::HECiphertext = 
    bootstrap(ct, scheme.oper, scheme.atk, scheme.rlk, scheme.boot)