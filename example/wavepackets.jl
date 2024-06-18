using GLMakie

# Definition of a wave packet (for the purposes of this example)
wp(ω, α, x0) = x -> cos(ω*x)*exp(-α*(x-x0)^2)

# A figure with sliders to show how the parameters affect the wave packet
pts = range(-10,10,2000)
vals(ω, α, x0) = wp(ω, α, x0).(pts)
f = Figure()
ax = Axis(f[1,1], limits=(-10,10,-2,2), title = L"f(x) = \cos(ω x) \cdot \exp(-α (x - x_0)^2)", titlesize=20)
sg = SliderGrid(f[2,1], 
               (label = "ω", range = range(0,100,101), format="{:.1f}", startvalue=10),
               (label = "α", range = 10 .^range(-2,2,101), format="{:.1f}", startvalue=1),
               (label = "x_0", range = range(-5,5,101), format="{:.1f}", startvalue=0),
               tellwidth=false)
parameters = [s.value for s ∈ sg.sliders]
lines!(pts, lift(vals, parameters...))

using ApproxFun, MTFun

# Coeffs plot figure

MT = MalmquistTakenaka()
f = Figure()
ax = Axis(f[1,1], limits=(0,1000,0,10), title = L"f(x) = \cos(ω x) \cdot \exp(-α (x - x_0)^2)", titlesize=20)
sg = SliderGrid(f[2,1], 
               (label = "ω", range = range(0,100,101), format="{:.1f}", startvalue=10),
               (label = "α", range = 10 .^range(-2,2,101), format="{:.1f}", startvalue=1),
               (label = "x_0", range = range(-5,5,101), format="{:.1f}", startvalue=0),
               tellwidth=false)
parameters = [s.value for s ∈ sg.sliders]
u = Fun(lift(wp,parameters...).val, MT)
lines!(abs.(coefficients(u)))

using Plots

Plots.plot(u - 1im*u)