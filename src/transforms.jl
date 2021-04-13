

## Transform
hasfasttransform(::MalmquistTakenaka) = true
points(S::MalmquistTakenaka,n::Int) = real(S.λ) .+ imag(S.λ)*tan.(range(0 ,length=n,step=π/n))
function plan_transform(S::MalmquistTakenaka,vals::AbstractVector)
    n = length(vals)
    weights = sqrt(π/abs(imag(S.λ))).*(points(S,n) .- conj(S.λ))
    TransformPlan(S,(weights,FastTransforms.plan_fft(complex(vals))),Val{false})
end
function plan_itransform(S::MalmquistTakenaka,cfs::AbstractVector)
    weights = sqrt(abs(imag(S.λ))/π) ./ (points(S,length(cfs)) .- conj(S.λ))
    ITransformPlan(S,(weights,FastTransforms.plan_ifft(complex(cfs))),Val{false})
end
function *(P::TransformPlan{T,S,false},vals::AbstractVector) where {T,S<:MalmquistTakenaka}
    weights,fftplan = P.plan
    _mt_reorder_scale_and_phase_shift(fftplan*(weights.*complex(vals)))
end
function *(P::ITransformPlan{T,S,false},cfs::AbstractVector) where {T,S<:MalmquistTakenaka}
    weights,ifftplan = P.plan
    weights.*(ifftplan*(_imt_reorder_and_phase_shift(complex(cfs))))
end

function _mt_reorder_scale_and_phase_shift(cfs::AbstractVector)
    n = length(cfs)
    tmp = FastTransforms.ifftshift(cfs)
    tmp .*= (one(eltype(cfs))*im).^(range(-ceil(n/2),length=n))/n
    newcfs = copy(cfs)
    newcfs[1:2:n] = tmp[div(n,2)+mod(n,2)+1:n]
    newcfs[2:2:n] = tmp[div(n,2)+mod(n,2):-1:1]
    newcfs
end

function _imt_reorder_and_phase_shift(cfs::AbstractVector)
    newcfs = copy(cfs)
    n = length(cfs)
    newcfs[div(n,2)+mod(n,2)+1:n] = cfs[1:2:n]
    newcfs[div(n,2)+mod(n,2):-1:1] = cfs[2:2:n]
    newcfs .*= (-one(eltype(cfs)*im)).^(range(-ceil(n/2),length=n))
    FastTransforms.fftshift(newcfs)
end
