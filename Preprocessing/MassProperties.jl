using LinearAlgebra

function totalmass(craft::basicaircraft)
    totalmass = 0.0
    for mass in craft.masses
        totalmass = totalmass + mass.mass
    end
    return totalmass
end

function centerofmass(craft::basicaircraft)
    centerofmass = zeros(3)
    for mass in craft.masses
        centerofmass = centerofmass + mass.mass*mass.pos
    end
    centerofmass = centerofmass*(1/totalmass(craft))
end

function skewsymetricsquared(pos)
    x = pos[1]
    y = pos[2]
    z = pos[3]
    D = zeros(3,3)
    D[1,1] = y^2 +z^2
    D[1,2] = -x*y
    D[1,3] = -x*z
    D[2,1] = -y*x
    D[2,2] = x^2 +z^2
    D[2,3] = -y*z
    D[3,1] = -z*x
    D[3,2] = -z*y
    D[3,3] = x^2+y^2
    return D
end

function rotateItocenter(mass::body)
    xdir = mass.dirx
    ydir = mass.diry
    zdir = cross(xdir,ydir)
    R = zeros(3,3)
    R[:,1] = xdir
    R[:,2] = ydir
    R[:,3] = zdir

    return R*mass.I*R'
end

function inertia(mass::pointmass,center)
    d = mass.pos - center
    J = mass.mass*skewsymetricsquared(d)
    return J
end

function inertia(mass::body,center)
    d = mass.pos-center
    D = skewsymetricsquared(d)
    J = rotateItocenter(mass)
    return J + mass.mass*D
end

function inertia(craft::basicaircraft)
     J = zeros(3,3)
     center = centerofmass(craft)
     for mass in craft.masses
        J = J + inertia(mass,center)
     end
     return J
end