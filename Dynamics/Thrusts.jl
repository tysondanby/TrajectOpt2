include("../Structs/Propulsion.jl")
using LinearAlgebra
function utothrust(propulsion::simplethrust,ui)
    T = propulsion.dir * propulsion.Tmax * ui
    return T
end

function utomoment(propulsion::simplethrust,ui)
    M = cross(propulsion.pos,utothrust(propulsion,ui))
    return M
end