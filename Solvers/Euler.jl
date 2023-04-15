function euler(f,x0,tspan)
    dt = 1e-3
    x = deepcopy(x0)
    xhist = []
    thist = []

    for t = tspan[1]:dt:tspan[2]
        push!(thist,t)
        push!(xhist,x)
        x = x + f(x,[],t)*dt
        #println(x)
    end
    return xhist,thist
end
