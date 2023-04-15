include("buildVLmodel.jl")
include("../MassProperties.jl")
using FLOWMath

function forcesVLMsteady(craft::basicaircraft, alp,bet)
    VLsurfs=buildVLmodel(craft)
    global testsurfs = VLsurfs
    Vinf = 1.0
    alpha = alp*pi/180 # angle of attack
    beta = bet*pi/180 # sideslip angle
    #println("AB:",alpha,beta)
    Omega = [0.0, 0.0, 0.0] # rotational velocity around the reference location
    #it is zero due to being steady
    fs = Freestream(Vinf, alpha, beta, Omega)

    Sref, cref = reference(craft) # reference area

    bref = Sref/cref # reference span
    rref = centerofmass(craft) # reference location for rotations/moments (typically the c.g.)
    Vinf = 1.0 # reference velocity (magnitude)
    ref = Reference(Sref, cref, bref, rref, Vinf)
    global testref = ref
    global testfs =fs
    system = steady_analysis(VLsurfs, ref, fs)
    global testsys = system
    CF, CM = body_forces(system; frame=Wind())
    CD, CY, CL = CF
    Cl, Cm, Cn = CM
    CD = far_field_drag(system)
    #NOTE: this is induced drag only
    return CD,CY,CL, Cl,Cm,Cn
end

function sweepforcesVLMsteady(craft::basicaircraft, alphas,betas)
    CDs = zeros(length(alphas),length(betas))
    CYs = similar(CDs)
    CLs = similar(CDs)
    Cls = similar(CDs)
    Cms = similar(CDs)
    Cns = similar(CDs)
    alphamax = 15.0
    alphamaxindex=indexin(maximum(filter(t -> t <alphamax,alphas)),alphas)[1]
    alphamin = -15.0 #TODO: find a way to avoid coding these in hard.
    alphaminindex=indexin(minimum(filter(t -> t >alphamin,alphas)),alphas)[1]
    betamax = 15.0   #      Possibly, I could use the wing polars.
    betamaxindex=indexin(maximum(filter(t -> t <betamax,betas)),betas)[1]
    betamin = -15.0
    betaminindex=indexin(minimum(filter(t -> t >betamin,betas)),betas)[1]
    
    for i = alphaminindex:1:alphamaxindex
        alphaused = alphas[i] 
        for j = betaminindex:1:betamaxindex
                betaused = betas[j]
                CDs[i,j],CYs[i,j],CLs[i,j], Cls[i,j],Cms[i,j],Cns[i,j] = forcesVLMsteady(craft,alphaused,betaused)
                #println(CDs[i,j],CYs[i,j],CLs[i,j], Cls[i,j],Cms[i,j],Cns[i,j])
        end
    end

    #top
    for i = alphamaxindex:1:length(alphas)
        CDs[i,:],CYs[i,:],CLs[i,:], Cls[i,:],Cms[i,:],Cns[i,:] = CDs[alphamaxindex,:],CYs[alphamaxindex,:],CLs[alphamaxindex,:], Cls[alphamaxindex,:],Cms[alphamaxindex,:],Cns[alphamaxindex,:]
    end
    #bottom
    for i = 1:1:alphaminindex
        CDs[i,:],CYs[i,:],CLs[i,:], Cls[i,:],Cms[i,:],Cns[i,:] = CDs[alphaminindex,:],CYs[alphaminindex,:],CLs[alphaminindex,:], Cls[alphaminindex,:],Cms[alphaminindex,:],Cns[alphaminindex,:]
    end
    for j = betamaxindex:1:length(betas)
        CDs[:,j],CYs[:,j],CLs[:,j], Cls[:,j],Cms[:,j],Cns[:,j] = CDs[:,betamaxindex],CYs[:,betamaxindex],CLs[:,betamaxindex], Cls[:,betamaxindex],Cms[:,betamaxindex],Cns[:,betamaxindex] 
    end
    for j = 1:1:betaminindex
        CDs[:,j],CYs[:,j],CLs[:,j], Cls[:,j],Cms[:,j],Cns[:,j] = CDs[:,betaminindex],CYs[:,betaminindex],CLs[:,betaminindex], Cls[:,betaminindex],Cms[:,betaminindex],Cns[:,betaminindex] 
    end

    return CDs,CYs,CLs, Cls,Cms,Cns
end
