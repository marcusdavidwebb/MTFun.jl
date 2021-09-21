function horner(cfs::AbstractVector,z::T) where T
    ret = zero(T)
    for k = length(cfs):-1:1
        ret = ret .* z + cfs[k]
    end
    ret
end

## Evaluate at a point using Horner's method (above)
function evaluate(cfs::AbstractVector,S::MalmquistTakenaka,x)
    T = eltype(x)
    z = (-one(T)im)*(S.λ .- x)./(conj(S.λ) .- x)
    ret = (-one(T)im)*horner(cfs[1:2:end],z) ./ (conj(S.λ) .- x)
    ret = ret .+ horner(cfs[2:2:end],inv.(z)) ./ (S.λ .- x)
    ret *= sqrt(abs(imag(S.λ))/π)
    ret
end

function evaluate(cfs::AbstractVector,S::Weideman,x)
    T = eltype(x)
    z = (-one(T)im)*(S.λ .- x)./(conj(S.λ) .- x)
    ret = horner(cfs[1:2:end],z)
    ret = ret .+ horner([zero(T);cfs[2:2:end]],inv.(z))
    ret
end

## Reorders coefficients to be indexed by ℤ
function hatify(cfs::AbstractVector)
    newcfs = copy(cfs)
    n = length(cfs)
    newcfs[div(n,2)+mod(n,2)+1:n] = cfs[1:2:n]
    newcfs[div(n,2)+mod(n,2):-1:1] = cfs[2:2:n]
    -div(n,2)-mod(n,2):div(n,2)-1, newcfs
end

@recipe function f(F::Fun{Sp,T,Vector{T}}) where {Sp<:Union{MalmquistTakenaka,Weideman},T}
    S = space(F)
    λ = S.λ
    a = real(λ) - 30*imag(λ)
    b = real(λ) + 30*imag(λ)
    cfs = coefficients(F)
    if length(cfs) < 200
        pad!(cfs,200)
    end
    n = length(cfs)
    pts = points(S,n)
    pts = [pts[div(n,2)+2:n]; pts[1:div(n,2)+mod(n,2)]]
    vals = values(F)
    vals = [vals[div(n,2)+2:n]; vals[1:div(n,2)+mod(n,2)]]
    
    inds = findall(x -> a ≤ x ≤ b,pts)
    @series begin 
        pts[inds],real(vals[inds])
    end
    if maximum(abs.(imag(vals[inds]))) ≥ 1e-6
        @series begin
            pts[inds],imag(vals[inds])
        end
    end
end

