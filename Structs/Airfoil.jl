abstract type Airfoil end

mutable struct NACA4 <: Airfoil
    number::Int16
    #polars are just vectors of equaly spaced coefficient values in [-pi:pi]
    clpolar::Vector{Any}
    cdpolar::Vector{Any}
    cmpolar::Vector{Any}
end

mutable struct NACA6 <: Airfoil
    number::Int32
    #polars are just vectors of equaly spaced coefficient values in [-pi:pi]
    clpolar::Vector{Any}
    cdpolar::Vector{Any}
    cmpolar::Vector{Any}
end

mutable struct splineairfoil <: Airfoil
    points::Vector{Any}  #X and Y points that trace out the airfoil
    #polars are just vectors of equaly spaced coefficient values in [-pi:pi]
    clpolar::Vector{Any}
    cdpolar::Vector{Any}
    cmpolar::Vector{Any}
end 