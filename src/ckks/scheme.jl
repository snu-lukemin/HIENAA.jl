"""
    CKKSScheme(sketch::CKKSParamSketch)
    CKKSScheme(param::CKKSParameters)
    CKKSScheme(oper::CKKSOperator)

CKKSScheme is a struct for the CKKS homomorphic encryption scheme. It contains the following fields: oper, entor, rlk, and atk.
- oper: CKKSOperator
- entor: Union{Missing,Encryptor}
- rlk: Union{Missing,RLEV}
- atk: Dict{Int64,RLEV}
"""
mutable struct CKKSScheme <: HEScheme
    oper::CKKSOperator
    boot::Union{Missing,CKKSBootstrapper}
    entor::Union{Missing,Encryptor}
    rlk::Union{Missing,RLEV}
    atk::Dict{Int64,RLEV}

    CKKSScheme(sketch::CKKSParamSketch)::CKKSScheme = CKKSScheme(CKKSParameters(sketch))

    function CKKSScheme(param::CKKSParameters)::CKKSScheme
        oper = CKKSOperator(param)

        new(oper, missing, missing, missing, Dict{Int64,RLEV}())
    end

    CKKSScheme(oper::CKKSOperator)::CKKSScheme = new(oper, missing, missing, missing, Dict{Int64,RLEV}())
end
