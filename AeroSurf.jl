include("Airfoil.jl")
abstract type AeroSurf end

mutable struct wing <: AeroSurf
    pos         #Position on Aircraft
    points::Vector{Any}      #Locations of each referenced airfoil. root to tip
    chorddist::Vector{Any}   #Vector of chord lengths
    twistdist::Vector{Any}
    airfoils::Vector{Airfoil}    #Vector of airfiol objects
    #TODO: write an inner constructor for creating a wing using just the dihedral and sweep instead of "points"
end

mutable struct bodyofrotation <: AeroSurf
    start       #Position of start of body
    dir::Vector{Any}         #Direction vector for centerline of body of rotation
    l::Vector{Any}           #vector of distances along dir
    r::Vector{Any}           #vector of radii corresponding to l
end
