mutable struct test
    a::Vector{Any}
    b::Int16
    c
end

mutable struct ODE <: Function 
    f
end
