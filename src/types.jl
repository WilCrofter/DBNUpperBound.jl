# custom data types 

""" NotBigReals

    A convenient type for trapping certain arguments. SpecialFunctions `gamma` and `zeta` can take real, but not complex, BigFloats.
    """
const NotBigInt = Union{Int64,Int32,Int16,Int8}
const NotBigFloat = Union{Float64, Float32, Float16}
const NotBigRational = Rational{T} where {T<:NotBigInt}
const NotBigReal = Union{NotBigInt, NotBigFloat, NotBigRational}

""" NotBigComplex

    A convenient type for trapping certain arguments. SpecialFunctions `lgamma`, `gamma` and `zeta` can take real, but not complex, BigFloats.
    """
const NotBigComplex = Complex{T} where {T <: NotBigReal}
