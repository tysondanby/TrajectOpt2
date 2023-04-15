include("Optimization/trajectmatch.jl")
include("buildCRC3.jl")
tspan = (0,3)
ts =  collect(tspan[1]:.5:tspan[2])
xs = similar(ts)
zs = similar(ts)
@. xs = 3*ts^2
@. zs = 10*ts^.5
x = LinearInterpolation((ts),xs)
z =LinearInterpolation((ts),zs)

function path(t)
    return [x(t),z(t)]
end
us0 = [] 
nupoints = 12
for i = 1:1:nupoints
    push!(us0,zeros(2))
end

pos0 = [0.0,0.0,0.0]
Theta0 = 90.0
v0 = [0.0,0.0,10.0]
x0 = [pos0...,  0.0,Theta0,0.0,   v0...,   0.0,0.0,0.0]

println("beginning optimization")

usstar, fopt = trajectmatch(path,tspan,us0,CRC3,x0)