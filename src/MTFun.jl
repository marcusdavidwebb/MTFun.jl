module MTFun

using ApproxFun, FastTransforms, Plots
import ApproxFun: points, plan_transform, plan_itransform, TransformPlan, ITransformPlan, domain, canonicaldomain, spacescompatible, evaluate, ConcreteMultiplication, ConcreteDerivative
import Base: *, first, last, getindex
export MalmquistTakenaka, hatify

"""
`MalmquistTakenaka(λ)` is the space spanned by
```
    (-i)^n sqrt(imag(λ)/π) * (λ - x)^n / (conj(λ) - x)^(n+1)
```
for all n ∈ ℤ.
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

end # module
