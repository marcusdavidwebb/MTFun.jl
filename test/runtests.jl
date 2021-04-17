using MTFun
using Test
using ApproxFun
using Plots

@testset "MTFun.jl" begin
    # Write your tests here.
    MT = MalmquistTakenaka()
    φ = (x,n) -> (1.0im)^n * sqrt(2/π) * (1+2im*x)^n / (1-2im*x)^(n+1)
    @test coefficients(Fun(x->φ(x,0),MT,12)) ≈ [1;zeros(11)] 
    @test coefficients(Fun(x->φ(x,5),MT,12)) ≈ [zeros(10);1;0]
    @test coefficients(Fun(x->φ(x,0),MT,11)) ≈ [1;zeros(4)] 
    @test coefficients(Fun(x->φ(x,5),MT,11)) ≈ [zeros(10);1]
    
    f = x -> 1/(1+x^2)

    @test f(pt) == φ0(pt)
end
