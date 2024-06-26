# Multiplication operators MT x MT to MT
function Multiplication(f::Fun{MT},sp::MT) where MT <: MalmquistTakenaka 
    if spacescompatible(space(f),sp)
        ConcreteMultiplication(f,sp)
    else
        error("Failure to build multiplication operator. Spaces not compatible.")
    end
end

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

    scaling = 1/sqrt(4*abs(imag(space(M).λ))*π)

    val = zero(eltype(coefficients(M.f)))
    if ind1 ≤ length(coefficients(M.f))
        val += coefficients(M.f)[ind1]*scaling
    end
    if ind2 ≤ length(coefficients(M.f))
        val -= 1im*coefficients(M.f)[ind2]*scaling
    end
    return val
end

function bandwidths(M::ConcreteMultiplication{MT,MT}) where MT <: MalmquistTakenaka
    n = length(coefficients(M.f))
    n = n + mod(n,2)
    return n,n
end

# Multiplication operators W x MT to MT and W x W to W (same matrix)
function Multiplication(f::Fun{W},sp::S) where {W<:Weideman, S <: Union{Weideman,MalmquistTakenaka}}  
    if spacescompatible(space(f),sp)
        ConcreteMultiplication(f,sp)
    else
        error("Failure to build multiplication operator. Spaces not compatible.")
    end
end

space(M::ConcreteMultiplication{W,S,T}) where {W<:Weideman, S <: Union{Weideman,MalmquistTakenaka}, T} = M.space
rangespace(M::ConcreteMultiplication{W,S,T}) where {W<:Weideman, S <: Union{Weideman,MalmquistTakenaka},T} = domainspace(M)

Base.stride(M::ConcreteMultiplication{U,V}) where {U<:Weideman,V<:Union{Weideman,MalmquistTakenaka}} =
    stride(M.f)

function getindex(M::ConcreteMultiplication{W,S},k::Integer,j::Integer) where {W<:Weideman, S <: Union{Weideman,MalmquistTakenaka}} 
    # convert k,j to biinfinite indices
    k=iseven(k) ? -k÷2 : (k-1)÷2
    j=iseven(j) ? -j÷2 : (j-1)÷2
    
    # convert k - j to original indices
    biinf_ind1 = k - j
    ind = (biinf_ind1 ≥ 0) ? 1 + 2biinf_ind1 : -2biinf_ind1

    val = zero(eltype(coefficients(M.f)))
    if ind ≤ length(coefficients(M.f)) 
        val += coefficients(M.f)[ind]
    end
    return val
end

function bandwidths(M::ConcreteMultiplication{W,S}) where {W <: Weideman, S <: Union{Weideman,MalmquistTakenaka}}
    n = length(coefficients(M.f))
    n = n + mod(n,2) - 1
    return n,n
end




## Derivative MT to MT

function Derivative(sp::MT,order::Integer) where MT <: MalmquistTakenaka
    ConcreteDerivative(sp,order)
end

space(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = D.space
rangespace(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = domainspace(D)
bandwidths(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = 2*D.order,2*D.order
Base.stride(D::ConcreteDerivative{MT}) where MT <: MalmquistTakenaka = 2*D.order

function getindex(D::ConcreteDerivative{MalmquistTakenaka{T},K,KK},k::Integer,j::Integer) where {T,K,KK}
    m = D.order
    bandind = j-k

    if m == 1
        if bandind == 0
            ((-one(T))^(j+1))im * (k-mod(k+1,2)) /(2imag(space(D).λ))
        elseif bandind == 2
            -div(k+1,2) /(2imag(space(D).λ))
        elseif bandind == -2
            div(j+1,2) /(2imag(space(D).λ))
        else
            zero(T)
        end
    elseif m == 2
        if bandind == 0
            -one(T) * (6div(k+1,2)^2 - 6div(k+1,2) + 2 ) /(2imag(space(D).λ))^2
        elseif bandind == 2
            ((-one(T))^k)*4im*div(k+1,2)^2 /(2imag(space(D).λ))^2
        elseif bandind == -2
            ((-one(T))^(j+1))*4im*div(j+1,2)^2 /(2imag(space(D).λ))^2
        elseif bandind == 4
            div(k+1,2)*div(k+3,2)/(2imag(space(D).λ))^2
        elseif bandind == -4
            div(j+1,2)*div(j+3,2)/(2imag(space(D).λ))^2
        else
            zero(T)
        end
    else
        error("Higher order Malmquist--Takenaka derivatives are not supported yet")
    end
end
