function horner(cfs::AbstractVector,z::T) where T
    ret = zero(T)
    for k = length(cfs):-1:1
        ret = ret .* z + cfs[k]
    end
    ret
end

## Evaluate at a point using Horner's method (above)
function evaluate(cfs::AbstractVector,S::MalmquistTakenaka,x::T) where T
    z = (-ones(T)im)*(S.λ .- x)./(conj(S.λ) .- x)
    ret = (-ones(T)im)*horner(cfs[1:2:end],z) ./ (conj(S.λ) .- x)
    ret += horner(cfs[2:2:end],inv.(z)) ./ (S.λ .- x)
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

function plot(f::Fun{MalmquistTakenaka{T},T,Vector{T}}) where T
    λ = space(f).λ
    a = real(λ) - 20*imag(λ)
    b = real(λ) + 20*imag(λ)
    plot([a,b],f)
end

function plot(interval,f::Fun{MalmquistTakenaka{T},T,Vector{T}}) where T
    cfs = coefficients(f)
    S = space(f)
    n = min(200,length(cfs))
    pts = points(S,n)
    inds = findall(x -> interval[1] ≤ x ≤ interval[2],pts)
    vals = itransform(S,cfs)
    plot(pts[inds],real(vals[inds]))
    if minimum(abs.(imag(vals))) ≥ 1e-4
        plot!(pts[inds],imag(vals[inds]))
    end
end

