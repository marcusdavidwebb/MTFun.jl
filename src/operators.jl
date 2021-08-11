# Multiplication

Base.stride(M::ConcreteMultiplication{U,V}) where {U<:MalmquistTakenaka,V<:MalmquistTakenaka} =
    stride(M.f)

function getindex(M::ConcreteMultiplication{MalmquistTakenaka{T},MalmquistTakenaka{T},T},k::Integer,j::Integer) where {T}
    tmp = (-1)^k * div(k,2) + (-1)^j * div(j,2)
    ind1 = 2*abs(tmp) + ((tmp >= 0) ? 1 : 0)
    tmp += 1
    ind2 = 2*abs(tmp) + ((tmp >= 0) ? 1 : 0)
    (coefficients(M.f)[ind1] - 1im*coefficients(M.f)[ind2]) / sqrt(2*convert(T,π))
 end


## Derivative

function Derivative(sp::MT,order::Integer) where MT <: MalmquistTakenaka
    if order == 1
        ConcreteDerivative(sp,order)
    else
        Derivative(sp,1)^order
    end
end

space(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = domainspace(D)
rangespace(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = domainspace(D)
bandwidths(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = 2*D.order,2*D.order
Base.stride(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = 2*D.order

function getindex(D::ConcreteDerivative{MalmquistTakenaka{T},K,KK},k::Integer,j::Integer) where {T,K,KK}
    if k == j
        ((-one(T))^(k+1))im * (k-mod(k+1,2)) /(2*imag(space(D).λ))
    elseif k == j+2
        #(-one(T))^(j+1) * 
        div(j+1,2) /(2*imag(space(D).λ))
    elseif j == k+2
        #(-one(T))^k * 
        -div(k+1,2) /(2*imag(space(D).λ))
    else
        zero(T)
    end
end
