# Multiplication

#Base.stride(M::ConcreteMultiplication{U,V}) where {U<:MalmquistTakenaka,V<:MalmquistTakenaka} =
#    stride(M.f)

#function getindex(M::ConcreteMultiplication{MalmquistTakenaka{T},MalmquistTakenaka{T},T},k::Integer,j::Integer) where {T}
#  tmp = (-1)^k * div(k,2) + (-1)^j * div(j,2)
#  ind1 = 2*abs(tmp) + ((tmp >= 0) ? 1 : 0)
#  tmp += 1
#  ind2 = 2*abs(tmp) + ((tmp >= 0) ? 1 : 0)
#  return (coefficients(M.f)[ind1] - 1im*coefficients(M.f)[ind2]) / sqrt(2*convert(T,π))
# end


## Derivative

#Derivative(sp::MalmquistTakenaka{T},order::Integer) where T =
#    ConcreteDerivative(sp,order)

#rangespace(D::ConcreteDerivative{MalmquistTakenaka{T}}) where T =
#    MalmquistTakenaka{T}()

#bandwidths(D::ConcreteDerivative{MalmquistTakenaka}) = -2*D.order,2*D.order
#Base.stride(D::ConcreteDerivative{MalmquistTakenaka}) = 2*D.order

#function getindex(D::ConcreteDerivative{MalmquistTakenaka{S},K,T},k::Integer,j::Integer) where {S,K,T}
    0
#end
