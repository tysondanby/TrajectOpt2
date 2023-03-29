abstract type Control end

struct simplecontroller <: Control
    wing #number corresponding to which wing surface the control surf is on
    dcl#/du
    dcd
    dcm
end