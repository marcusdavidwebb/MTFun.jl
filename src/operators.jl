# Multiplication
Multiplication(f::Fun{MT},sp::MT) where MT <: MalmquistTakenaka = ConcreteMultiplication(f,sp)

space(M::ConcreteMultiplication{MT,MT,T}) where {MT <: MalmquistTakenaka, T} = M.space
rangespace(M::ConcreteMultiplication{MT,MT,T}) where {MT <: MalmquistTakenaka,T} = domainspace(M)

Base.stride(M::ConcreteMultiplication{U,V}) where {U<:MalmquistTakenaka,V<:MalmquistTakenaka} =
    stride(M.f)

function getindex(M::ConcreteMultiplication{MT,MT},k::Integer,j::Integer) where {MT <: MalmquistTakenaka} 
    # convert k,j to biinfinite indices
    k=iseven(k) ? -k÷2 : (k-1)÷2
    j=iseven(j) ? -j÷2 : (j-1)÷2
    
    # convert k - j and k - j - 1 to original indices
    biinf_ind1 = k - j
    ind1 = (biinf_ind1 ≥ 0) ? 1 + 2biinf_ind1 : -2biinf_ind1
    biinf_ind2 = k - j - 1
    ind2 = (biinf_ind2 ≥ 0) ? 1 + 2biinf_ind2 : -2biinf_ind2

    val = zero(eltype(coefficients(M.f)))
    if ind1 ≤ length(coefficients(M.f)) 
        val += coefficients(M.f)[ind1] / sqrt(2π)
    end
    if ind2 ≤ length(coefficients(M.f))
        val -= 1im*coefficients(M.f)[ind2] / sqrt(2π)
    end
    return val
end

function bandwidths(M::ConcreteMultiplication{MT,MT}) where MT <: MalmquistTakenaka
    n = length(coefficients(M.f))
    n = n + mod(n,2)
    return n,n
end

## Derivative

function Derivative(sp::MT,order::Integer) where MT <: MalmquistTakenaka
    if order == 1
        ConcreteDerivative(sp,order)
    else
        Derivative(sp,1)^order
    end
end

space(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = D.space
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
