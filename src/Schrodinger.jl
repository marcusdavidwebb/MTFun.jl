
function plan_R0_step(u::Fun,V,ε,Δt)
    itransplan = plan_itransform(space(u),coefficients(u))
    Vmultiplier = exp.(-1im*(Δt/ε)*V.(points(u)))
    transplan = plan_transform(space(u),values(u))
    return itransplan, Vmultiplier, transplan
end

function R0_step(u::Fun, itransplan, Vmultiplier, transplan)
    vals = itransplan*coefficients(u)
    newcoeffs = transplan*(Vmultiplier.*vals)
    return Fun(space(u),newcoeffs)
end

function plan_R1_step(u::Fun,ε,Δt,order = 0)
    n = length(coefficients(u))
    D2 = Derivative(space(u),2)[1:n,1:n]
    if order == 2
        expD2top = 2I + 1im*ε*Δt*D2
        expD2bottom = 2I - 1im*ε*Δt*D2
    elseif order == 4
        expD2top = 12I + 6im*ε*Δt*D2 - (ε*Δt*D2)^2
        expD2bottom = 12I - 6im*ε*Δt*D2 - (ε*Δt*D2)^2
    elseif order == 6
        expD2top = 120I + 60im*ε*Δt*D2 - 12*(ε*Δt*D2)^2 - 1im*(ε*Δt*D2)^3
        expD2bottom = 120I - 60im*ε*Δt*D2 - 12*(ε*Δt*D2)^2 + 1im*(ε*Δt*D2)^3
    elseif order == 0
        if bandwidths(D2) == (0,0) # for example, Fourier
            expD2top = complex(D2)
            expD2top.data = exp.(1im*ε*Δt*expD2top.data)
        else
            expD2top = exp(Matrix(1im*ε*Δt*D2))
        end
        expD2bottom = I
    end
    return (expD2bottom, expD2top)
end

function R1_step(u::Fun,expD2)
    newcoeffs = expD2[1] \ (expD2[2] * coefficients(u))
    return Fun(space(u),newcoeffs)
end

function plan_R2_step(u0::Fun,V1,V2,ε,Δt)
    n = length(coefficients(u0))
    D2term = (Δt^3 * ε /6 ) * Derivative(space(u0),2)[1:n,1:n]
    itransplan = plan_itransform(space(u0),coefficients(u0))
    transplan = plan_transform(space(u0),values(u0))
    V1term = (Δt^3 / 12ε) * V1.(points(u0)).^2
    V2term = V2.(points(u0))
    return D2term, itransplan, transplan, V1term, V2term
end

function R2_step(u::Fun, D2term, itransplan, transplan, V1term, V2term, Krylov_dim=6)
    q1 = complex(coefficients(u))
    T = eltype(q1)
    normq1 = norm(q1)
    q1 = q1/normq1
    n = length(q1)
    m = Krylov_dim
    Q = zeros(T,n,m) 
    Q[:,1] = q1
    H = zeros(T,m,n)
        
    for k=1:m
        z = one(T)im * D2term * (transplan * (V2term .* (itransplan * Q[:,k])))
        z += one(T)im * transplan * (V2term .* (itransplan * (D2term * Q[:,k])))
        z +=  one(T)im * transplan * (V1term .* (itransplan * Q[:,k]))
        for i=1:k
            H[i,k] = Q[:,i]'*z
            z = z - H[i,k]*Q[:,i]
        end
        if k < m
            H[k+1,k] = norm(z)
            if H[k+1,k] == zero(T)
                return
            end
            Q[:,k+1] = z / H[k+1,k]
        end
    end
    return Q, H[1:m,1:m], normq1*(Q*exp(H[1:m,1:m]))[:,1]
end

function plan_Lie_Trotter_step(u::Fun,V,ε,Δt,order=0)
    itransplan, Vmultiplier, transplan = plan_R0_step(u,V,ε,Δt)
    expD2 = plan_R1_step(u,ε,Δt,order)
    return itransplan, Vmultiplier, transplan, expD2
end

function Lie_Trotter_step(u::Fun, itransplan, Vmultiplier, transplan, expD2)
    newu = R0_step(u,itransplan, Vmultiplier, transplan)
    return R1_step(newu,expD2)
end

function Lie_Trotter_evolve(u0::Fun,V,ε,Δt,T,order=0)
    M = Integer(floor(T/Δt))
    u0 = Fun(space(u0),complex(coefficients(u0)))
    retFuns = Vector{typeof(u0)}(undef,M+1)
    retFuns[1] = u0
    itransplan, Vmultiplier, transplan, expD2 = plan_Lie_Trotter_step(u,V,ε,Δt,order)
    for k = 1:M
        retFuns[k+1] = Lie_Trotter_step(retFuns[k], itransplan, Vmultiplier, transplan, expD2)
    end
    return 0:Δt:T, retFuns
end


function plan_Strang_step(u::Fun,V,ε,Δt,order=0)
    itransplan, Vmultiplier, transplan = plan_R0_step(u,V,ε,Δt/2)
    expD2 = plan_R1_step(u,ε,Δt,order)
    return itransplan, Vmultiplier, transplan, expD2
end

function Strang_step(u::Fun, itransplan, Vmultiplier, transplan, expD2)
    newu = R0_step(u,itransplan, Vmultiplier, transplan)
    newu = R1_step(newu,expD2)
    return R0_step(newu,itransplan, Vmultiplier, transplan)
end

function Strang_evolve(u0::Fun,V,ε,Δt,T,order=0)
    M = Integer(floor(T/Δt))
    u0 = Fun(space(u0),complex(coefficients(u0)))
    retFuns = Vector{typeof(u0)}(undef,M+1)
    retFuns[1] = u0
    itransplan, Vmultiplier, transplan, expD2 = plan_Strang_step(u0,V,ε,Δt,order)
    for k = 1:M
        retFuns[k+1] = Strang_step(retFuns[k], itransplan, Vmultiplier, transplan, expD2)
    end
    return 0:Δt:T, retFuns
end





  
struct HermiteFSE{T<:Real} <: Space{Line{false,T},T}
    t :: T
end

HermiteFSE{T}() where T = HermiteFSE{T}(zero(T))
HermiteFSE() = HermiteFSE{Float64}()

domain(::HermiteFSE{T}) where T = Line{false,T}()
canonicaldomain(::HermiteFSE{T}) where T = Line{false,T}()
normalization(::Type{T}, ::HermiteFSE, ::Int) where T = one(T)

spacescompatible(a::HermiteFSE,b::HermiteFSE) = a.t == b.t
canonicalspace(::HermiteFSE{T}) where T = HermiteFSE(zero(T))

function Derivative(sp::HermiteFSE,order::Integer)
    ConcreteDerivative(sp,order)
end

space(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = D.space
rangespace(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = domainspace(D)
bandwidths(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = D.order,D.order
Base.stride(D::ConcreteDerivative{HFSE}) where HFSE <: HermiteFSE = D.order

function getindex(D::ConcreteDerivative{HermiteFSE{T},K,KK},k::Integer,j::Integer) where {T,K,KK}
    m = D.order
    bandind = j-k

    if m == 1
        if bandind == 1
            convert(T,sqrt(k/2))
        elseif bandind == -1
            convert(T,-sqrt(j/2))
        else
            zero(T)
        end
    elseif m == 2
        if bandind == 0
            convert(T,.5-k)
        elseif bandind == 2
            convert(T,sqrt(k*(k+1))/2)
        elseif bandind == -2
            convert(T,sqrt(j*(j+1))/2)
        else
            zero(T)
        end
    else
        error("Higher order HermiteFSE derivatives are not supported yet")
    end
end

function evaluate(cfs::AbstractVector,::HermiteFSE,x)
    T = eltype(x)
    n = length(cfs)
    ret = zero(T)
    if (n > 0) ret += cfs[1] * exp(-x^2 / 2) * convert(T,π)^(-1/4) end
    if (n > 1) ret += cfs[2] * x * exp(-x^2 / 2) * (4/convert(T,π))^(1/4) end
    
    if (n > 2)
        hkm1 = convert(T,π)^(-1/4)       #h_0(x)
        hk = sqrt(one(T)*2) * x * hkm1  #h_1(x)
        sum_log_scale = zero(T)
        for k = 2:n-1
            # hk = h_k(x), hkm1 = h_{k-1}(x) (actually, recaled versions thereof)
            hkm1, hk = hk, sqrt(one(T)*2/k)*x*hk - sqrt((k-one(T))/k)*hkm1
            # rescale values
            scale = (x->(x<one(T)) ? one(T) : inv(x))(abs(hk))
            hk *= scale
            hkm1 *= scale
            # keep track of final rescaling factor
            sum_log_scale += log(scale)
            ret += cfs[k+1] * hk * exp(-x^2 / 2 - sum_log_scale)
        end
    end
    return ret
end

function eval_Hermite_function(n::Integer,x::T) where T <: Number
    # Stabalised version of 3 term recurrence (see appendix of B. Bunck, BIT Numerical Mathematics (2009))
    # See also https://www.numbercrunch.de/blog/2014/08/calculating-the-hermite-functions/
        
    if (n == 0)
        return exp(-x^2 / 2) * convert(T,π)^(-1/4)
    elseif (n == 1)
        return x * exp(-x^2 / 2) * (4/convert(T,π))^(1/4)
    else
        hkm1 = convert(T,π)^(-1/4)       #h_0(x)
        hk = sqrt(one(T)*2) * x * hkm2  #h_1(x)
        sum_log_scale = zero(T)
        for k = 2:n 
            # hk = h_k(x), hkm1 = h_{k-1}(x), hkm2 = h_{k-2}(x) (actually, recaled versions thereof)
            hkm1, hk = hk, sqrt(one(T)*2/k)*x*hk - sqrt((k-one(T))/k)*hkm1
            # rescale values
            scale = (x->(x<one(T)) ? one(T) : inv(x))(abs(hk))
            hk *= scale
            hkm1 *= scale
            # keep track of final rescaling factor
            sum_log_scale += log(scale)
        end
        return hk * exp(-x^2 / 2 - sum_log_scale)
    end
end

## Transforms for HermiteFSE space
hasfasttransform(::MalmquistTakenaka) = true
# Gauss-Hermite points
points(::HermiteFSE,n::Int) = gausshermite(n)[1]

function slow_plan_Hermite_transform(n::Integer)

    slowC = zeros(n+1,n+1)
    x,~ = gausshermite(n+1)

    ψn = sqrt(n+1)*eval_Hermite_function.(n,x)

    for k = 0:n
        slowC[k+1,:] = eval_Hermite_function.(k,x)./ψn
    end

    valweights = 1 ./ ψn

    return valweights, slowC
end

function plan_Hermite_transform(n::Integer)

    C = zeros(n+1,n+1)
    x,~ = gausshermite(n+1)

    hkm1 = ones(n+1)*π^(-1/4)       #h_0(x)
    hk = sqrt(2) * x .* hkm1  #h_1(x)
    cum_log_scale = zeros(n+1,n+1)

    C[1,:] = hkm1
    if (n ≥ 1) C[2,:] = hk end

    for k = 2:n 
        # hk = h_k(x), hkm1 = h_{k-1}(x) (actually, recaled versions thereof)
        # reassign values
        hkm1, hk = hk, sqrt(2/k)* x .* hk - sqrt((k-1)/k)*hkm1
        # rescale values
        scale = (x->(x<1) ? 1 : inv(x)).(abs.(hk))
        hk .*= scale
        hkm1 .*=  scale
        C[k+1,:] = hk
        # keep track of rescaling factors
        cum_log_scale[k+1,:] = cum_log_scale[k,:] + log.(scale)
    end
    
    # the function values should be weighted by these
    valweights = exp.(cum_log_scale[n+1,:] + x.^2 / 2) ./ (sqrt(n+1)*C[n+1,:])
    
    # C contains rescaled Hermite polynomials. Need to convert to unscaled ψ_k(x_j)/(sqrt(n+1)ψ_n(x_j)).
    for k = 1:n+1
        C[k,:] = C[k,:]./C[n+1,:] .* exp.(cum_log_scale[n+1,:]-cum_log_scale[k,:] .- .5*log(n+1))
    end

    return valweights, C
end

function plan_transform(S::HermiteFSE,vals::AbstractVector)
    valweights, C = plan_Hermite_transform(length(vals)-1)
    TransformPlan(S,(valweights,C),Val{false})
end
function plan_itransform(S::HermiteFSE,cfs::AbstractVector)
    valweights, C = plan_Hermite_transform(length(cfs)-1)
    ITransformPlan(S,(valweights,C),Val{false})
end

function *(P::TransformPlan{T,S,false},vals::AbstractVector) where {T,S<:HermiteFSE}
    valweights,C = P.plan
    C*(valweights.*vals)
end
function *(P::ITransformPlan{T,S,false},cfs::AbstractVector) where {T,S<:HermiteFSE}
    valweights,C = P.plan
    (C' * cfs) ./ valweights
end