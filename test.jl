using Plots
include("Structs/Aircraft.jl")
include("Preprocessing/Polars.jl")
alphas = -180.0:0.5:180.0
res = 100000:100000:200000
emptypolar = zeros(length(alphas),length(res))
defaultairfoil = NACA4(2412,emptypolar,emptypolar,emptypolar,0.0,0.0)
origin = [0.0,0.0,0.0]
defaultchord = [1.0,1.0]
defaulttwist = [0.0,0.0]
airfoils = [defaultairfoil,defaultairfoil]
defaultwing = wing(origin,[origin,[0.0,1.0,0.0]],defaultchord,defaulttwist,airfoils,0.0)

wingpolars!(defaultwing,res,alphas)

plot(alphas,defaultwing.airfoils[1].clpolar[:,1])