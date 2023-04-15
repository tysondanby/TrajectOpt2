include("Structs/Aircraft.jl")
include("Structs/Atmosphere.jl")
include("Dynamics/Simulate.jl")

function dummy(x)
    return 0
end
#SET UP AIRCRAFT 
xle =  0.0
alphas = -180.0:0.5:180.0
res = 100000:100000:200000
emptypolar = zeros(length(alphas),length(res))
defaultairfoil = NACA4("0012",emptypolar,emptypolar,emptypolar,0.0,0.0)
origin = [0.0,0.0,0.0]
defaultchord = [0.08604,0.08604]
defaulttwist = [0.0,0.0]
airfoils = [defaultairfoil,defaultairfoil]
topwing = wing([xle,0,0.125],[[0,-2.952, 0],[0,2.952, 0]],defaultchord,defaulttwist,airfoils,0.50798)#TODO: S isn't calculated automatically.
bottomwing = wing([xle,0,-0.125],[[0,-2.952, 0],[0,2.952, 0]],defaultchord,defaulttwist,airfoils,0.50798)#TODO: S isn't calculated automatically.

surfs = [topwing, bottomwing]

#mass1 = pointmass([.1,0,0],.5)
#mass2 = pointmass([0,0.1,0],.5)
#mass3 = pointmass([0,0,0.1],.5)
#mass4 = pointmass([0.0,0,0],.5)
J = [
    100.0 0.0 0.0;
    0.0 0.0111 0.0;
    0.0 0.0 100.0
]
masses = [body([.02,0.0,0.0],[1.0,0,0],[0,1.0,0],1.36,J)]#[mass1,mass2,mass3,mass4]

pos = [0.0,0.0,0.125]
topthrust = simplethrust(-pos,[1.0,0.0,0.0],15.0,"Top")
bottomthrust = simplethrust(pos,[1.0,0.0,0.0],15.0,"Bottom")
prop = [topthrust, bottomthrust]

control = []#Not implemented in current version

statenames =  ["x","y","z","phi","theta","psi","u","v","w","p","q","r"]#x = x y z phi theta psi u v w p q r
controlnames = ["no control!"]
paramnames = ["none"]
dynamics = DAEdynamics(statenames,paramnames,controlnames,dummy)

S = 0.50798
c = 0.08604

standardATM = StationaryUniformAtmosphere(1.225)
CRC3 = basicaircraft(surfs,masses,prop,control,dynamics,S,c,"Rock",false)
CRC3.dynamics.f = basicdynamics2DODE(CRC3,standardATM)

println("CRC3 Built")

#SET UP SIMULATION
#Topu = [.5, .5, .5]
#Bottomu = [.5, .5, .5]
us = [[0., 1.],[0., 1.],[0.0,0.0]]
pos0 = [0.0,0.0,0.0]
Theta0 = 90.0
v0 = [0.0,0.0,10.0]
x0 = [pos0...,  0.0,Theta0,0.0,   v0...,   0.0,0.0,0.0]#x = x y z phi theta psi u v w p q r
dx0= [0.0,0.0,0.0,  0.0,0.0,0.0,   0.0,0.0,0.0,   0.0,0.0,0.0]#filler now
tspan = (0.0,5.0)
println("preprocessing")
xhist,uhist,thist = simulate(CRC3,us,dx0,x0,tspan)#us is a vector of u-vectors