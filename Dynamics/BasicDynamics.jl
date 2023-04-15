include("../Preprocessing/VortexLattice/tabulateVLM.jl")
include("Thrusts.jl")
include("../Structs/Atmosphere.jl")
include("ReferenceFrames.jl")

using Interpolations

function basicdynamics2DODE(craft,atm)
    mass = totalmass(craft)
    forces, moments = getforcesmoments(craft,atm)
    J = inertia(craft)
    return function f(x,u,params,t)
        xdot = zeros(length(x))
        #x = x y z phi theta psi u v w p q r
        #phi = x[4]
        #theta = x[5]
        #psi = x[6]
        
        #omega = [phi, theta, psi]
        #p = x[10]
        q = x[11]
        #r = x[12]
        #println("Forces",forces(xdot,x,u,params,t))
        #println("Moments",moments(xdot,x,u,params,t))

        xdot[1:3] = [x[7],x[8],x[9]]
        xdot[4] = 0.0#p + (q*sind(phi) + r*cosd(phi))*tand(theta)
        xdot[5] = q*180/pi#*cosd(phi) - r*sind(phi)
        xdot[6] = 0.0#(q*sind(phi) + r*cosd(phi))/cosd(theta)
        xdot[7:9]= forces(xdot,x,u,params,t)/mass  #x,y,z
        xdot[11] = moments(xdot,x,u,params,t)[2]/J[2,2]

        #2D
        xdot[8] = 0.0
        return xdot
    end
end

function basicdynamics3DODE(craft,atm)
    mass = totalmass(craft)
    forces, moments = getforcesmoments(craft,atm)
    J = inertia(craft)
    return function f(x,u,params,t)
        xdot = zeros(length(x))
        #x = x y z phi theta psi u v w p q r
        phi = x[4]
        theta = x[5]
        psi = x[6]
        
        omega = [phi, theta, psi]
        p = x[10]
        q = x[11]
        r = x[12]
        println("Forces",forces(xdot,x,u,params,t))
        println("Moments",moments(xdot,x,u,params,t))

        xdot[1:3] = [x[7],x[8],x[9]]
        xdot[4] = p + (q*sind(phi) + r*cosd(phi))*tand(theta)
        xdot[5] = q*cosd(phi) - r*sind(phi)
        xdot[6] = (q*sind(phi) + r*cosd(phi))/cosd(theta)
        xdot[7:9]= forces(xdot,x,u,params,t)/mass  #x,y,z
        xdot[10:12] = J \ (moments(xdot,x,u,params,t) - cross(omega, J*omega))
        return xdot
    end
end

function basicdynamics3D(craft,atm)
    mass = totalmass(craft)
    forces, moments = getforcesmoments(craft,atm)
    return function f(xdot,x,u,params,t)
        R = zeros(length(x))
        #x = x y z phi theta psi u v w p q r
        vdot = [xdot[7],xdot[8],xdot[9]]
        phi = x[4]
        theta = x[5]
        psi = x[6]
        phidot = xdot[10]
        thetadot = xdot[11]
        psidot = xdot[12]
        omegadot = [xdot[10],xdot[11],xdot[12]]
        omega = [phi, theta, psi]
        p = x[10]
        q = x[11]
        r = x[12]



        #Accelerations
        R[1:3]= forces(xdot,x,u,params,t) - mass*vdot  #x,y,z
        R[4:6]= I \ (moments(xdot,x,u,params,t) - cross(omega, I*omega)) - omegadot

        #Velocities
        R[7:9]= [xdot[1],xdot[2],xdot[3]] - [x[7],x[8],x[9]]
        R[10] = p + (q*sin(phi) + r*cos(phi))*tan(theta) -phidot
        R[11] = q*cos(phi) - r*sin(phi) -thetadot
        R[12] = (q*sin(phi) + r*cos(phi))/cos(theta) -psidot
        return R
    end
end

function getforcesmoments(craft::basicaircraft,atm::StationaryUniformAtmosphere)
    forcestabulated = false #TODO: make this check if CDs,CYs,CLs, Cls,Cms,Cns exist
    alphas = -90:0.5:90
    betas = -180:0.5:180
    CDs = zeros(length(alphas),length(betas))
    CYs = similar(CDs)
    CLs = similar(CDs)
    Cls = similar(CDs)
    Cms = similar(CDs)
    Cns = similar(CDs)
    

    if forcestabulated == false#TODO:else fill in from tabulated data.
        CDs,CYs,CLs, Cls,Cms,Cns = sweepforcesVLMsteady(craft, alphas,betas)
    end

    CDfunc = LinearInterpolation((alphas,betas),CDs)
    CLfunc = LinearInterpolation((alphas,betas),CLs)
    CYfunc = LinearInterpolation((alphas,betas),CYs)
    Clfunc = LinearInterpolation((alphas,betas),Cls)
    Cmfunc = LinearInterpolation((alphas,betas),Cms)
    Cnfunc = LinearInterpolation((alphas,betas),Cns)
    #global testCD = CDfunc
    #global testCL = CLfunc
    #global testCY = CYfunc
    #global testCl = Clfunc
    #global testCm = Cmfunc
    #global testCn = Cnfunc

    #global testCDs = CDs
    #global testCLs = CLs
    #global testCYs = CYs
    #global testCls = Cls
    #global testCms = Cms
    #global testCns = Cns

    COM=centerofmass(craft)
    S = craft.S
    c = craft.c
    Rho = atm.rho
    #Local functions
    
    function thrust(u,x)
        T = zeros(3)
        if length(u) != 0
            for i = 1:1:length(craft.propulsion)
                #println(utothrust(craft.propulsion[i],u[i]))
                T = T + utothrust(craft.propulsion[i],u[i])
            end
        end
        return bodytoinertial(T,x[4],x[5],x[6])
    end
    function gravity(m)
        return [0,0,-9.8]*m
    end
    function aeroforces(x)
        vinf = [-x[7],-x[8],-x[9]] #+ wind(x,atm)
        bodyx = bodytoinertial([1.0,0,0],x[4],x[5],x[6])
        bodyy = bodytoinertial([0,1.0,0],x[4],x[5],x[6])
        #bodyz = bodytoinertial([0,0,1.0],x[10],x[11],x[12])

        #Define drag,lift, Fy directions
        zwind = normalize(cross(vinf,bodyy))
        ywind= normalize(cross(zwind,vinf))
        liftdir = normalize(cross(vinf,bodyy))
        fydir = normalize(cross(liftdir,vinf))#TODO: make sure is right direction.  might be oposite.

        CDparasitic = 0.0151 #TODO: update this
        beta = atand(dot(normalize(vinf),bodyy),dot(ywind,bodyy))
        alpha = acosd(dot(zwind,bodyx)) - 90.0
        #println("AlphaBeta:",alpha,beta)
        #push!(AB,[alpha, beta])

        qinf = 0.5 * Rho * dot(vinf,vinf)
        CL = CLfunc(alpha,beta)
        CD = CDfunc(alpha,beta) + CDparasitic
        CY = CYfunc(alpha,beta)
        L = S*qinf*CL
        D = S*qinf*CD
        Fy = S*qinf*CY#TODO: make sure this works before 3Ding it up.

        Force = L * -liftdir + D * normalize(vinf) + Fy* fydir#Vinf of aircraft

        return Force#bodytoinertial(Force,x[4],x[5],x[6])
    end
    function thrustmoments(u,x)
        M = zeros(3)
        if length(u) != 0
            for i = 1:1:length(craft.propulsion)
                M = M + cross(craft.propulsion[i].pos-COM,utothrust(craft.propulsion[i],u[i]))
            end
        end
        return bodytoinertial(M,x[4],x[5],x[6])
    end
    function aeromoments(x)
        vinf = [-x[7],-x[8],-x[9]] #+ wind(x,atm)
        bodyx = bodytoinertial([1.0,0,0],x[4],x[5],x[6])
        bodyy = bodytoinertial([0,1.0,0],x[4],x[5],x[6])
        #bodyz = bodytoinertial([0,0,1.0],x[10],x[11],x[12])

        #Define wind directions
        zwind = normalize(cross(vinf,bodyy))
        ywind= normalize(cross(zwind,vinf))#TODO: make sure is right direction.  might be oposite.

        beta = atand(dot(normalize(vinf),bodyy),dot(ywind,bodyy))
        alpha = acosd(dot(zwind,bodyx)) - 90.0
        #println("AlphaBeta:",alpha,beta)

        qinf = 0.5 * Rho * dot(vinf,vinf)
        Cl = Clfunc(alpha,beta)
        Cm = Cmfunc(alpha,beta) #TODO: check my assumption: positive nose up
        Cn = Cnfunc(alpha,beta)
        l = S*qinf*Cl*c#TODO: verify if these are the right reference distances (i.e. is c the right distance?)
        m = S*qinf*Cm*c
        n = S*qinf*Cn*c

        moment = l*normalize(vinf) + m*ywind + n*zwind#Vinf of aircraft
        
        return moment#bodytoinertial(moment,x[4],x[5],x[6])
    end

    #Returned functions
    function forces(xdot,x,u,p,t)
        #println("Thrust:",thrust(u,x))
        #println("Aero Force:",aeroforces(x))
        return aeroforces(x) + gravity(totalmass(craft)) + thrust(u,x)
    end
    function moments(xdot,x,u,p,t)
        #println("Thrust Moment:",thrustmoments(u,x))
        #println("Aero Moment:",aeromoments(x))
        return aeromoments(x) + thrustmoments(u,x)
    end
    #global testF = forces
    #global testM = moments
    return forces,moments
end