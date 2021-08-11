module MTFun

using ApproxFun, ApproxFunBase, FastTransforms, Plots, LinearAlgebra
import ApproxFun: points, plan_transform, plan_itransform, TransformPlan, ITransformPlan, domain, canonicaldomain, spacescompatible, evaluate
import ApproxFunBase: ConcreteMultiplication, Derivative, ConcreteDerivative, space, rangespace, bandwidths
import Base: *, first, last, getindex
export MalmquistTakenaka, hatify

"""
`MalmquistTakenaka(λ)` is the space spanned by
```
    φ_n(x) = (-i)^{n+1} sqrt(|imag(λ)|/π) (λ - x)^n / (λ̄ - x)^(n+1)
```
for all n ∈ ℤ. The canonical case λ = i/2 can also be written
```
    φ_n(x) =  i^n sqrt(2/π) (1 + 2ix)^n / (1- 2ix)^(n+1)
``` 
"""

struct MalmquistTakenaka{T<:Complex} <: Space{Line{false,T},T}
    λ :: T
end

MalmquistTakenaka{T}() where T = MalmquistTakenaka{T}(one(T)*im/2)
MalmquistTakenaka() = MalmquistTakenaka{ComplexF64}()

domain(::MalmquistTakenaka{T}) where T = Line{false,T}()
canonicaldomain(::MalmquistTakenaka{T}) where T = Line{false,T}()
normalization(::Type{T}, ::MalmquistTakenaka, ::Int) where T = one(T)

Base.first(::Fun{<:MalmquistTakenaka{T}}) where T = zero(T)
Base.last(::Fun{<:MalmquistTakenaka{T}}) where T= zero(T)

spacescompatible(a::MalmquistTakenaka,b::MalmquistTakenaka) = a.λ == b.λ
canonicalspace(S::MalmquistTakenaka) = S

include("transforms.jl")
include("evalandplot.jl")
include("operators.jl")
include("Schrodinger.jl")

end # module
