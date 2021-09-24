module MTFun

using ApproxFun, FastTransforms, Plots, LinearAlgebra, FastGaussQuadrature
import ApproxFun: points, plan_transform, plan_itransform, TransformPlan, ITransformPlan, domain, canonicaldomain, spacescompatible, evaluate
import ApproxFun: Multiplication, ConcreteMultiplication, Derivative, ConcreteDerivative, space, rangespace, bandwidths
import Base: *, first, last, getindex
export MalmquistTakenaka, Weideman, HermiteFSE
export hatify

"""
`MalmquistTakenaka(λ)` is the space spanned by
```
    φₙ(x) = (-i)^{n+1} sqrt(|imag(λ)|/π) (λ - x)^n / (λ̄ - x)^(n+1)
```
for all n ∈ ℤ. The canonical case λ = i/2 can also be written
```
    φₙ(x) =  i^n sqrt(2/π) (1 + 2ix)^n / (1 - 2ix)^(n+1)
``` 
This basis is orthonormal in L₂(ℝ).
"""

struct MalmquistTakenaka{T<:Complex} <: Space{Line{false,T},T}
    λ :: T
end

MalmquistTakenaka{T}() where T = MalmquistTakenaka{T}(one(T)*im/2)
MalmquistTakenaka() = MalmquistTakenaka{ComplexF64}()

domain(::MalmquistTakenaka{T}) where T = Line{false,T}()
canonicaldomain(::MalmquistTakenaka{T}) where T = Line{false,T}()
normalization(::Type{T}, ::MalmquistTakenaka, ::Int) where T = one(T)

spacescompatible(a::MalmquistTakenaka,b::MalmquistTakenaka) = a.λ == b.λ
canonicalspace(S::MalmquistTakenaka) = S

"""
`Weideman(λ)` is the space spanned by
```
    σₙ(x) = (-i)^n (λ - x)^n / (λ̄ - x)^n
```
for all n ∈ ℤ. The canonical case λ = i/2 can also be written
```
    σₙ(x) =  i^n (1 + 2ix)^n / (1 - 2ix)^n
``` 
This basis is orthonormal in L₂(ℝ,w), where w(x) = |imag(λ)| / (π*(x^2 - 2Real(λ)x + |λ|^2))
"""

struct Weideman{T<:Complex} <: Space{Line{false,T},T}
    λ :: T
end

Weideman{T}() where T = Weideman{T}(one(T)*im/2)
Weideman() = Weideman{ComplexF64}()

domain(::Weideman{T}) where T = Line{false,T}()
canonicaldomain(::Weideman{T}) where T = Line{false,T}()
normalization(::Type{T}, ::Weideman, ::Int) where T = one(T)

spacescompatible(a::Weideman,b::Weideman) = a.λ == b.λ
canonicalspace(S::Weideman) = S

spacescompatible(a::Weideman,b::MalmquistTakenaka) = a.λ == b.λ
spacescompatible(a::MalmquistTakenaka,b::Weideman) = a.λ == b.λ

include("transforms.jl")
include("evalandplot.jl")
include("operators.jl")
include("Schrodinger.jl")

end # module
