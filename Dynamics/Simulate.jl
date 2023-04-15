using DifferentialEquations
using NLsolve
include("BasicDynamics.jl")
include("../Solvers/DiscreteMethods.jl")
include("../Solvers/Euler.jl")

function simulate(craft,us,dx0,x0,tspan)
    dt = (tspan[2]-tspan[1])/(length(us)-1)
    ts = collect(tspan[1]:dt:tspan[2])
    interpu = LinearInterpolation(ts,us)
    if length(interpu) == 0
        interpu = (x)-> [0.0]
    end
    #TODO: only currently works with DAEs
    dynamics = injectus(craft.dynamics.f,interpu,tspan)
    #println("Preprocessing complete. Solving...")#TEMP
    #dx0 = firststep(dynamics,x0,tspan[1])
    #global testdx0 = dx0
    #global testinterpu = interpu
    #global testdynam = dynamics
    #dt = 0.001 #TODO: make this determined somewhere
    #tsteps = collect(tspan[1]:dt:tspan[2])
    #prob = ODEProblem(dynamics,x0,tspan,reltol = 1e-6)
    #xhist = discreteEuler(dynamics,[],x0,tsteps)#TODO: use something better than euler
    #prob = DAEProblem(dynamics, dx0, x0, tspan)#, differential_vars = differential_vars)
    #global testprob = prob
    xhist,thist = euler(dynamics,x0,tspan)
    #sol = solve(prob)
    #xhist = sol.u
    #thist = sol.t
    uhist = @. interpu(thist)
    return mapreduce(permutedims,vcat,xhist),uhist,thist
end

function injectusDAE(f,interpu,tspan)#takes a function in my form, inputs the us, and returns an equation in diffeq form
    #dt = (tspan[2]-tspan[1])/(length(us)-1)
    #ts = collect(tspan[1]:dt:tspan[2])
    #interpu = LinearInterpolation(ts,us)
    return function diffeqform(out,xdot,x,p,t)
        out[:] = f(xdot,x,interpu(t),p,t)
        return out
    end
end

function injectus(f,interpu,tspan)#takes a function in my form, inputs the us, and returns an equation in diffeq form
    #dt = (tspan[2]-tspan[1])/(length(us)-1)
    #ts = collect(tspan[1]:dt:tspan[2])
    #interpu = LinearInterpolation(ts,us)
    return function diffeqform(x,p,t)
        xdot = zeros(length(x))
        xdot[:] = f(x,interpu(t),p,t)
        return xdot
    end
end

function firststep(f,x0,t)
    resids = deepcopy(x0)
    function R!(resids,dx0)
        resids=f(resids,dx0,x0,[],t)
        return resids
    end
    #solve resids = R(dx0) for dx0 that makes resids zero.
    results = nlsolve(R!,dx0)
    return results.zero
end