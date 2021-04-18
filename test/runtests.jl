using MTFun
using Test
using ApproxFun

@testset "Constructors, Evaluation, Transforms" begin
    λ = -0.7+0.9im
    MT = MalmquistTakenaka(λ)
    φ = (x,n) -> (-1.0im)^(n+1)*sqrt(abs(imag(λ))/π) * (λ - x)^n / (conj(λ) - x)^(n+1)
    @test coefficients(Fun(x->φ(x,0),MT,12)) ≈ [1;zeros(11)] 
    @test coefficients(Fun(x->φ(x,5),MT,12)) ≈ [zeros(10);1;0]
    @test coefficients(Fun(x->φ(x,0),MT,11)) ≈ [1;zeros(10)] 
    @test coefficients(Fun(x->φ(x,5),MT,11)) ≈ [zeros(10);1]
    
    f = x -> (1.0-2im)/(1+x^2)
    F = Fun(f,MT)
    x = 0.78
    @test F(x) ≈ f(x)

    f = x -> exp(-10x^2)
    F = Fun(f,MT)
    x = 0.78
    @test real(F(x)) ≈ f(x)

    vals = φ.(points(MT,10),0)
    cfs = [1.0+0.0im;zeros(9)]
    trans = MTFun.plan_transform(MT,vals)
    itrans = MTFun.plan_itransform(MT,cfs)
    @test cfs ≈ trans*vals
    @test vals ≈ itrans*cfs

    v = [0.83+0.0im;0.78;-1.13;0.65;1.11]
    trans = MTFun.plan_transform(MT,v)
    itrans = MTFun.plan_itransform(MT,v)
    @test v ≈ itrans*(trans*v)
    @test v ≈ trans*(itrans*v)
end