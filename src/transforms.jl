

## Transforms for MalmquistTakenaka Space
hasfasttransform(::MalmquistTakenaka) = true
# note that one of these points is always infinity
points(S::MalmquistTakenaka{T},n::Int) where T = real(S.λ) .+ imag(S.λ)*tan.(range(zero(real(T)),length=n,step=π/n))

function plan_transform(S::MalmquistTakenaka,vals::AbstractVector)
    weights = sqrt(π/abs(imag(S.λ))).*(conj(S.λ) .- points(S,length(vals)))
    TransformPlan(S,(weights,FastTransforms.plan_fft(complex(vals))),Val{false})
end
function plan_itransform(S::MalmquistTakenaka,cfs::AbstractVector)
    weights = sqrt(abs(imag(S.λ))/π) ./ (conj(S.λ) .- points(S,length(cfs)))
    ITransformPlan(S,(weights,FastTransforms.plan_bfft(complex(cfs))),Val{false})
end
function *(P::TransformPlan{T,S,false},vals::AbstractVector) where {T,S<:MalmquistTakenaka}
    weights,fftplan = P.plan
    _mt_reorder_scale_and_phase_shift(fftplan*(weights.*complex(vals)))
end
function *(P::ITransformPlan{T,S,false},cfs::AbstractVector) where {T,S<:MalmquistTakenaka}
    weights,ifftplan = P.plan
    vals = weights.*(ifftplan*(_imt_reorder_and_phase_shift(complex(cfs))))
    if (iseven(length(vals))) vals[div(length(vals),2)+1] = zero(eltype(vals)) end
    return vals
end

## Transforms for Weideman Space
hasfasttransform(::Weideman) = true
# note that one of these points is always infinity
points(S::Weideman{T},n::Int) where T = real(S.λ) .+ imag(S.λ)*tan.(range(zero(real(T)) ,length=n,step=π/n))

function plan_transform(S::Weideman,vals::AbstractVector)
    TransformPlan(S,FastTransforms.plan_fft(complex(vals)),Val{false})
end
function plan_itransform(S::Weideman,cfs::AbstractVector)
    ITransformPlan(S,FastTransforms.plan_bfft(complex(cfs)),Val{false})
end
function *(P::TransformPlan{T,S,false},vals::AbstractVector) where {T,S<:Weideman}
    fftplan = P.plan
    _mt_reorder_scale_and_phase_shift(fftplan*((-one(T)im)*complex(vals)))
end
function *(P::ITransformPlan{T,S,false},cfs::AbstractVector) where {T,S<:Weideman}
    ifftplan = P.plan
    (one(T)im)*ifftplan*(_imt_reorder_and_phase_shift(complex(cfs)))
end


# Deals with the fact that the positively and negatively indexed 
# coefficients need to be interlaced
function _mt_reorder_scale_and_phase_shift(cfs::Vector{T}) where T
    n = length(cfs)
    tmp = FastTransforms.fftshift(cfs)
    tmp .*= (-one(T)*im).^(range(-floor(n/2)-one(T),length=n))/n
    newcfs = copy(cfs)
    newcfs[1:2:n] = tmp[div(n,2)+1:n]
    newcfs[2:2:n] = tmp[div(n,2):-1:1]
    newcfs
end

# Deals with the fact that the positively and negatively indexed 
# coefficients need to be interlaced
function _imt_reorder_and_phase_shift(cfs::Vector{T}) where T
    n = length(cfs)
    newcfs = copy(cfs)
    newcfs[div(n,2)+1:n] = cfs[1:2:n]
    newcfs[div(n,2):-1:1] = cfs[2:2:n]
    newcfs .*= (one(T)*im).^(range(-floor(n/2)-one(T),length=n))
    FastTransforms.ifftshift(newcfs)
end
