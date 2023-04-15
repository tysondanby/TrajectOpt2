using Rotations
function bodytoinertial(v,phi,theta,psi)
    rot =  [phi,theta,psi]*pi/180
    vnew = RotZ(-rot[3])*v
    vnew = RotY(-rot[2])*vnew
    vnew = RotX(-rot[1])*vnew
    return vnew
end