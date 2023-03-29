abstract type Propulsion end

#T = Tmax*u
struct simplethrust <: Propulsion
    pos
    dir
    Tmax
    inputname #name of u
end


struct propthrust <: Propulsion
    pos
    dir
    inputname
    propeller
    motor
end

