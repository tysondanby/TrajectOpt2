abstract type Mass end

mutable struct pointmass <: Mass
    pos::Vector{Any}
    mass::Any
end

mutable struct body <: Mass
    pos::Vector{Any}
    dirx::Vector{Any} #Direction of body x axis relative to the Aircraft
    diry::Vector{Any} #Direction of body y axis relative to the Aircraft
    mass::Any
    I#::Array{Any,Any}
end