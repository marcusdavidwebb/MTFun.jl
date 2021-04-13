function horner(cfs::AbstractVector{T},z::T) where T
    ret = zero(T)
    for k = length(cfs):-1:1
        ret = ret .* z + cfs[k]
    end
    ret
end

## Evaluate at a point using Horner's method (above)
function evaluate(cfs::AbstractVector{T},S::MalmquistTakenaka{T},x::T) where T
    z = -one(T)*im*(S.λ .- x)./(conj(S.λ) .- x)
    ret = horner(cfs[1:2:end],z) ./ (conj(S.λ) .- x)
    ret -= 1.0im*horner(cfs[2:2:end],inv.(z)) ./ (S.λ .+ x)
    ret *= sqrt(abs(imag(S.λ))/π)
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
