include("AeroSurf.jl")
include("Mass.jl")
include("Propulsion.jl")
include("Control.jl")
include("Dynamics.jl")

abstract type Aircraft end

mutable struct basicaircraft <: Aircraft
    surfs::Vector{AeroSurf}       #vector of AeroSurf
    masses::Vector{Mass}      #vector of mass
    propulsion::Vector{Propulsion}  #Vector of propulsion objects
    control::Vector{Control}
    dynamics::Dynamics    #A function  for the dynamics
    S #Reference area
    c
    name
    reanalyzeaero
end