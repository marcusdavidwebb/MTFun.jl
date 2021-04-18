
Complex{ComplexF64}
isa(ComplexF64,)
using ApproxFun
Line{false,ComplexF64}

ComplexF64 <: Complex

using FastTransforms
v = complex(rand(BigFloat,10))
plan_fft(v)

eltype(1.)

isa(v,AbstractVector{Complex{BigFloat}})

using Plots, LaTeXStrings

pts = LinRange(-10,10,1000)
MT = MalmquistTakenaka()

f1 = x -> 1/(1+x^2)
plot(pts,f1.(pts),legend=false,linewidth=2)
  annotate!(6,.7,text(L"\frac{1}{1+x^2}",30))
savefig("f1.pdf")

f2 = x-> 1/(1+x^2)^3
plot(pts,Inf.*pts,legend=false)
  plot!(pts,f2.(pts),linewidth=2)
  annotate!(6,.7,text(L"\frac{1}{(1+x^2)^3}",30))
savefig("f2.pdf")

f3 = x -> cos(x)/(1+x^2)
plot(pts,Inf.*pts,legend=false)
  plot!(pts,Inf.*pts)
  plot!(pts,f3.(pts),linewidth=2)
  annotate!(6,.7,text(L"\frac{\cos(x)}{1+x^2}",30))
savefig("f3.pdf")

cols = Plots.palette(:default)

g1 = Fun(f1,MT,1000)
inds1,cfs1 = hatify(coefficients(g1))
plot(inds1,abs.(cfs1),yscale=:log10,legend = false,xlims=(-70,69))
annotate!(35,1e-3,text(L"\color[rgb]{0,0.6056,0.97}\mathcal{O}'\left(3^{-|n|}\right)",20),ylims=(1e-20,1000))
savefig("g1.pdf")


g2 = Fun(f2,MT,1000)
inds2,cfs2 = hatify(coefficients(g2))
plot(inds1,abs.(cfs1),yscale=:log10,legend = false)
  plot!(inds2,abs.(cfs2),yscale=:log10,xlims=(-70,69),ylims=(1e-20,1000))
annotate!(50,1e-7,text(L"\color[rgb]{0,0.6056,0.97}\mathcal{O}'\left(3^{-|n|}\right)",15))
  annotate!(50,1e-3,text(L"\color[rgb]{0.8888,0.43565,0.278123}\mathcal{O}'\left(3^{-|n|}\right)",15))
savefig("g2.pdf")

g3 = Fun(f3,MT)
inds3,cfs3 = hatify(coefficients(g3))
plot(inds1,abs.(cfs1),yscale=:log10,legend = false)
  plot!(inds2,abs.(cfs2),yscale=:log10)
  plot!(inds3,abs.(cfs3),yscale=:log10,xlims=[-250,249],ylims=(1e-20,1000))
annotate!(400,1e-12,text(L"\color[rgb]{0,0.6056,0.97}\mathcal{O}'\left(3^{-|n|}\right)",15))
  annotate!(400,1e-8,text(L"\color[rgb]{0.8888,0.43565,0.278123}\mathcal{O}'\left(3^{-|n|}\right)",15))
  annotate!(400,.25e-1,text(L"\color[rgb]{0.24222,0.643275,0.30445}\mathcal{O}\left(|n|^{-5/4}\right)",15))
savefig("g3.pdf")
savefig("MTcoefficients.png")

f4 = x -> exp(-x^2)
plot(pts,f4.(pts),legend=false,linewidth=2)
  annotate!(6,.7,text(L"\mathrm{e}^{-x^2}",30))
savefig("f4.pdf")

f5 = x-> exp(-(x-3)^2)
plot(pts,f5.(pts),color=cols[2],legend=false,linewidth=2)
  annotate!(-3,.7,text(L"\mathrm{e}^{-(x-3)^2}",30))
savefig("f5.pdf")

f6 = x -> cos(20x)*exp(-x^2)
plot(pts,f6.(pts),color=cols[3],legend=false,linewidth=2)
  annotate!(6,.7,text(L"\cos(20 x) \, \mathrm{e}^{-x^2}",22))
savefig("f6.pdf")


g4 = Fun(f4,MT,2000)
inds4,cfs4 = hatify(coefficients(g4))
plot(inds4,abs.(cfs4),yscale=:log10,legend = false,xlims=(-500,499))
  annotate!(300,1e-3,text(L"\color[rgb]{0,0.6056,0.97}\mathcal{O}'\left(\rho^{-|n|^{2/3}}\right)",20),ylims=(1e-20,1000))
savefig("g4.pdf")

g5 = Fun(f5,MT,3000)
inds5,cfs5 = hatify(coefficients(g5))
plot(inds4,abs.(cfs4),yscale=:log10,legend = false)
  plot!(inds5,abs.(cfs5),yscale=:log10,xlims=(-1000,999),ylims=(1e-20,1000))
  annotate!(650,1e-7,text(L"\mathcal{O}'\left(\rho^{-|n|^{2/3}}\right)",15))
savefig("g5.pdf")

g6 = Fun(f6,MT,4000)
inds6,cfs6 = hatify(coefficients(g6))
plot(inds4,abs.(cfs4),yscale=:log10,legend = false)
  plot!(inds5,abs.(cfs5),yscale=:log10)
  plot!(inds6,abs.(cfs6),yscale=:log10,xlims=[-1500,1499],ylims=(1e-20,1000))
  annotate!(1100,1e-7,text(L"\mathcal{O}'\left(\rho^{-|n|^{2/3}}\right)",15))
savefig("g6.pdf")








g4 = Fun(f4,MT,2000)
  inds4,cfs4 = hatify(coefficients(g4))
  plot(inds4,max.(eps(),abs.(cfs4)),yscale=:log10,legend = false,xlims=[-500,499])

g5 = Fun(f5,MT,4000)
  inds5,cfs5 = hatify(coefficients(g5))
  plot(inds4,max.(eps(),abs.(cfs4)),yscale=:log10,legend = false)
    plot!(inds5,max.(eps(),abs.(cfs5)),yscale=:log10,xlims=[-500,499])

g6 = Fun(f6,MT)
  inds6,cfs6 = hatify(coefficients(g6))
  plot(inds4,max.(eps(),abs.(cfs4)),yscale=:log10,legend = false)
    plot!(inds5,max.(eps(),abs.(cfs5)),yscale=:log10)
    plot!(inds6,max.(eps(),abs.(cfs6)),yscale=:log10,xlims=[-500,499])








cfs = [0.5;0.5im]
g = Fun(S,cfs)
plot(pts,real(g.(pts)))
plot!(pts,imag(g.(pts)))
f = x -> 1/(1+4x^2) * sqrt(2/pi)
plot!(pts,real(f.(pts)))
  plot!(pts,imag(f.(pts)))

g = Fun(S,[1.0;0.0])
plot(pts,real(g.(pts)))
  plot!(pts,imag(g.(pts)))
f = x -> sqrt(2/pi) / (1-2im*x)
plot!(pts,real(f.(pts)))
  plot!(pts,imag(f.(pts)))

g = Fun(S,[0.0;1.0])
  plot(pts,real(g.(pts)))
    plot!(pts,imag(g.(pts)))
f = x -> (-1.0im)*sqrt(2/pi) / (1+2im*x)
  plot!(pts,real(f.(pts)))
    plot!(pts,imag(f.(pts)))



f = x -> 1/(1+x^2)
F = Fun(f,MalmquistTakenaka())
S = MalmquistTakenaka()
Multiplication(F::Fun{MalmquistTakenaka{T}},S::MalmquistTakenaka{T}) where T = ConcreteMultiplication{MalmquistTakenaka{T},MalmquistTakenaka{T},T}(F,S)

M1 = Multiplication(F,S)
typeof(M1)
@which M1[1,1]
@which getindex(M1,3,4)

G = Fun(f,Chebyshev())

M2 = Multiplication(G,Chebyshev())

@which M[1,1]

differentiate(F)
ApproxFun.Derivative(MalmquistTakenaka())
canonicaldomain(MalmquistTakenaka())
domain(MalmquistTakenaka())
canonicalspace(MalmquistTakenaka())

conversion_type(canonicalspace(MalmquistTakenaka()),sp)

plot(F)

g = x-> exp(x)
@which Fun(g)
H = differentiate(G)


@which ApproxFun.DefaultDerivative(Chebyshev())



plot(G)

using MTFun, ApproxFun, Plots
MT = MalmquistTakenaka()
φ = (x,n) -> (1.0im)^n * sqrt(2/π) * (1+2im*x)^n / (1-2im*x)^(n+1)
coefficients(Fun(x->φ(x,0),MT,12))
coefficients(Fun(x->φ(x,-1),MT,12))
coefficients(Fun(x->φ(x,1),MT,12))
coefficients(Fun(x->φ(x,-2),MT,12))
coefficients(Fun(x->φ(x,2),MT,12))
coefficients(Fun(x->φ(x,-3),MT,12))
coefficients(Fun(x->φ(x,3),MT,12))
coefficients(Fun(x->φ(x,0),MT,3))
coefficients(Fun(x->φ(x,-1),MT,5))
coefficients(Fun(x->φ(x,1),MT,5))
coefficients(Fun(x->φ(x,-2),MT,5))
coefficients(Fun(x->φ(x,2),MT,5))

f = x-> 1/(1+x^2)

F = Fun(f,MT,1000)
F(0.3+0.0im)
f(0.3)
F(0.3)
plot(F)

F = Fun(f,MT)
plot(hatify(coefficients(F))[1],abs.(hatify(coefficients(F))[2]),yscale=:log10)

f = x-> exp(-x^2)

F = Fun(f,MT)



vals = φ.(points(MT,10),0)
cfs = [1.0+0.0im;zeros(9)]
trans = MTFun.plan_transform(MT,vals)
itrans = MTFun.plan_itransform(MT,cfs)

cfs - trans*vals

vals - itrans*cfs

λ = randn()+randn()im
MT = MalmquistTakenaka(λ)
φ = (x,n) -> (-1.0im)^(n+1)*sqrt(abs(imag(λ))/π) * (λ - x)^n / (conj(λ) - x)^(n+1)
coefficients(Fun(x->φ(x,0),MT,12))
coefficients(Fun(x->φ(x,-1),MT,12))
coefficients(Fun(x->φ(x,1),MT,12))
coefficients(Fun(x->φ(x,-2),MT,12))
coefficients(Fun(x->φ(x,2),MT,12))
coefficients(Fun(x->φ(x,-3),MT,12))
coefficients(Fun(x->φ(x,3),MT,12))
coefficients(Fun(x->φ(x,0),MT,3))
coefficients(Fun(x->φ(x,-1),MT,5))
coefficients(Fun(x->φ(x,1),MT,5))
coefficients(Fun(x->φ(x,-2),MT,5))
coefficients(Fun(x->φ(x,2),MT,5))

f = x->1/(1+x^4)

F = Fun(f,MT)
length(coefficients(F))
plot(hatify(coefficients(F))[1],abs.(hatify(coefficients(F))[2]),yscale=:log10)


