using Ipopt, SNOW
include("../Dynamics/Simulate.jl")

function trajectmatch(path,tspan,us0,craft,x0)#finds the set of us that best match a given path
    designVars, nu = ustodesignvars(us0)
    function objective!(g,designvars)
        us = designvarstous(designvars,nu)
        g = 0.0*g
        return trajectmatchobj(us,path,tspan,craft,x0)
    end
    
    ng = 1
    lx = -ones(length(designVars))  # lower bounds on x
    ux = ones(length(designVars))  # upper bounds on x
    lg = -ones(ng)
    ug = ones(ng)
    ip_options = Dict(
        "tol" => 1e-3,
        "max_iter" => 300
        )
    solver = IPOPT(ip_options)
    options = Options(derivatives =CentralFD();solver)#TODO: for some reason the optimizer evaluates the objective as zero unless finite diff is used.
    #global dv = designVars
    #global testobjective = objective!
    xout, fopt, info = minimize(objective!, designVars, ng, lx, ux, lg, ug,options)
    usstar = designvarstous(xout,nu)
    return usstar, fopt
end


function trajectmatchobj(us,path,tspan,craft,x0)
    dx0 = zeros(length(x0))
    xhist,uhist,thist = simulate(craft,us,dx0,x0,tspan)
    ts = Vector{Float64}(thist)
    xs = xhist[:,1]
    zs = xhist[:,3]
    #deltaxs = similar(xs)
    #deltazs = similar(zs)
    error = similar(zs)
    global testx = xs
    global testz = zs
    global testt = thist
    global testpath = path
    for i = 1:1:length(thist)
        deltaxs = xs[i] - path(ts[i])[1]
        deltazs = zs[i] - path(ts[i])[2]
        error[i] = deltaxs^2 + deltazs^2
    end
    obj = trapz(ts,error)/(tspan[2]-tspan[1])
    return obj
end

function trapz(X::T, Y::T) where {T <: AbstractVector}
    # Check matching vector length
    @assert length(X) == length(Y)
    # Initialize Output
    out = similar(X)
    out[1] = 0
    # Iterate over arrays
    for i = 2:1:length(X)
      out[i] = out[i-1] + 0.5*(X[i] - X[i-1])*(Y[i] + Y[i-1])
    end
    # Return output
    return out[end]
end

function ustodesignvars(us)
    designvars = Vector{typeof(us[1][1])}()
    nu=length(us[1])
    for i = 1:1:length(us)
        for j = 1:1:nu
            push!(designvars,us[i][j])
        end
    end
    return designvars, nu
end

function designvarstous(designvars,nu)
    us = Vector{Vector{typeof(designvars[1])}}()
    steps = trunc(Int,length(designvars)/nu)
    for i = 1:1:steps
        newu = Vector{typeof(designvars[1])}()
        for j = 1:1:nu
            push!(newu,designvars[nu*(i-1)+j])
        end
        push!(us,newu)
    end
    return us
end
