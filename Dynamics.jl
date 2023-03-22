abstract type Dynamics end

mutable struct ODEdynamics <: Dynamics # xdot = f(x,u,p,t) 
    statenames::Vector{String}
    controlnames::Vector{String}
    paramnames::Vector{String}
    f::Function
end

mutable struct DAEdynamics <: Dynamics # 0 = f(xdot,x,u,p,t) 
    statenames::Vector{String}
    controlnames::Vector{String}
    paramnames::Vector{String}
    f::Function
end
