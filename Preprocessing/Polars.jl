include("../Structs/Aircraft.jl")
include("Geometry.jl")
using Xfoil, LinearAlgebra
#using Plots #TODO: TEMP

#Generates polars for each airfoil on each wing. Finds Sref for each wing
function aircraftpolars!(craft::basicaircraft,res,alphas)
    for i = 1:1:length(craft.surfs)
        if typeof(craft.surfs[i])== wing
            wingpolars!(craft.surfs[i],res,alphas)
        end
    end
end

function wingpolars!(w::wing,res,alphas)
    #Compute AR
    S = 0.0
    for i =2:1:length(w.points)
        S = S + abs(w.points[i][2]-w.points[i-1][2])*0.5*(w.chorddist[i]+w.chorddist[i-1])
    end
    S = S *2 #TODO: only works for symetric wings
    b = abs(norm(w.points[end]-w.points[1]))*2 #TODO: only works for symetric wings
    AR = b^2/S
    #Assign polars to each section of the wing
    for i = 1:1:length(w.airfoils)
        airfoilpolars!(w.airfoils[i],res,alphas,AR)
    end
    w.Sref = S/2
end

function airfoilpolars!(foil::NACA4,res,alphas,AR)
    #TODO First look up if this has already been computed
    x,y = NACA4points(foil.number)
    foil.clpolar, foil.cdpolar, foil.cmpolar, foil.clmax, foil.clmin = pointstopolars(x,y,res,alphas,AR)
end

function airfoilpolars!(foil::NACA6,res,alphas,AR)
    x,y = NACA6points(foil.number)
    foil.clpolar, foil.cdpolar, foil.cmpolar, foil.clmax, foil.clmin = pointstopolars(x,y,res,alphas,AR)
end

function airfoilpolars!(foil::splineairfoil,res,alphas,AR)
    x = []
    y = []
    for i = 1:1:length(foil.points)
        push!(x,foil.points[i][1])
        push!(y,foil.points[i][2])
    end
    foil.clpolar, foil.cdpolar, foil.cmpolar, foil.clmax, foil.clmin = pointstopolars(x,y,res,alphas,AR)
end

#This is the entire polar from -180 to 180
#Returns 3 2D polars (alpha and Re)
function pointstopolars(x,y,res,alphas,AR)
    clmax = 0.0
    clmin = 0.0
    da = alphas[2]-alphas[1]
    na = length(alphas)
    Xfoil.set_coordinates(x,y)
    xr, yr = Xfoil.pane()
    
    #alpha = -180:0.1:180 # range of angle of attacks, in degrees
    #re = 1e5 # Reynolds number
    clpolar = zeros(na,length(res))
    cdpolar = zeros(na,length(res))
    cmpolar = zeros(na,length(res))


    startsweep = 0
    endsweep = 0
    for i = 1:1:na
        if alphas[i] < -35
            startsweep = i
        elseif alphas[i] < 35
            endsweep = i
        end
    end
    #startsweep = 1
    #endsweep = na


    for i = 1:1:length(res)
        c_l = zeros(na)
        c_d = zeros(na)
        c_m = zeros(na)
        converged =zeros(Bool,na)
        #println(xr)
        #println(yr)
        c_l[startsweep:endsweep], c_d[startsweep:endsweep], c_dp, c_m[startsweep:endsweep], converged[startsweep:endsweep] = Xfoil.alpha_sweep(xr, yr, alphas[startsweep:endsweep], res[i], iter=100, zeroinit=false, printdata=false, reinit=true)
        c_l, c_d,  c_m, clmax, clmin = fillplate!(c_l, c_d, c_m, converged,alphas,AR)#fillvit!(c_l, c_d, c_m, converged,alphas,AR)
        #Assign results for this Re to the polars.
        clpolar[:,i] = c_l
        cdpolar[:,i] = c_d
        cmpolar[:,i] = c_m
    end

    return clpolar, cdpolar, cmpolar, clmax, clmin
end

#Fills in viterna calulation for a polar. returns 3 1D polars (alphas varied at constant Re)
function fillvit!(cls,cds,cms,converged,alphas,AR)
    clmax = 0.0
    alphaclmax = 0.0
    clmin = 0.0
    alphaclmin = 0.0
    cdmax = 0.0
    cdclmax = 0.0
    cdclmin = 0.0
    imin = 1
    imax = 1
    for i = 1:1:length(alphas)
        if converged[i] == true
            if cls[i] > clmax
                clmax = cls[i]
                alphaclmax = alphas[i]
                cdclmax = cds[i]
                imax = i
            elseif cls[i] < clmin
                clmin = cls[i]
                alphaclmin = alphas[i]
                cdclmin = cds[i]
                imin = i
            end
        end
    end
    #Calc viturna
    if AR >= 50
        cdmax = 2.01
    else
        cdmax = 1.11 + 0.018*AR
    end
    Ap = (clmax - cdmax*sind(alphaclmax)*cosd(alphaclmax))*sind(alphaclmax)/(cosd(alphaclmax)^2)
    Bp = cdclmax - cdmax*sind(alphaclmax)^2/cosd(alphaclmax)
    Am = (-clmin - cdmax*sind(-alphaclmin)*cosd(-alphaclmin))*sind(-alphaclmin)/(cosd(-alphaclmin)^2)
    Bm = cdclmin - cdmax*sind(-alphaclmin)^2/cosd(-alphaclmin)
    for i = 1:1:length(alphas)
        if converged[i] == false && alphas[i] > 0.0
            newcl = 0.5*cdmax*sind(2*alphas[i])+Ap*(cosd(alphas[i])^2)/sind(alphas[i])
            if newcl < clmin || newcl*0 == NaN
                cls[i]=2*pi(alphas[i]-180)*pi/180
            else
                cls[i] =newcl
            end
            cds[i] = cdmax *sind(alphas[i])^2 + Bp*cosd(alphas[i])
            cms[i] = cms[imax] #TODO: Probably not quite right, but shouldn't be an issue if stall not reached.
        elseif converged[i] == false && alphas[i] < 0.0
            newcl = -1*(0.5*cdmax*sind(2*-alphas[i])+Ap*(cosd(-alphas[i])^2)/sind(-alphas[i]))
            if newcl > -clmin || newcl*0 == NaN
                cls[i]=2*pi(alphas[i]+180)*pi/180
            else
                cls[i] = newcl
            end
            cds[i] = cdmax *sind(-alphas[i])^2 + Bp*cosd(-alphas[i])
            cms[i] = cms[imin] #TODO: Probably not quite right, but shouldn't be an issue if stall not reached.
        end
    end
    

    return cls,cds,cms
end

function fillplate!(cls,cds,cms,converged,alphas,AR)
    clmax = 0.0
    clmin = 0.0
    imin = 1
    imax = 1
    for i = 1:1:length(alphas)
        if converged[i] == 1
            if cls[i] > clmax
                clmax = cls[i]
                imax = i
            elseif cls[i] < clmin
                clmin = cls[i]
                imin = i
            end
        end
    end
    convergedend = 0
    i = trunc(Int,length(alphas)/2)
    #println(converged)
    failed = 0
    #i = i +10
    while failed == 0
        if converged[i] == false
            convergedend = i-1
            failed = 1
        end
        i = i +1
    end
    convergedbegin = 0
    i = trunc(Int,length(alphas)/2)
    failed = 0
    #i = i -10
    while failed == 0
        if converged[i] == false
            convergedbegin = i+1
            failed = 1
        end
        i = i -1
    end

    for i = 1:1:length(alphas)
        if i > convergedend || i<convergedbegin
            cls[i] =2*sind(alphas[i])*cosd(alphas[i])
            cds[i] = 2*sind(alphas[i])^2
            if alphas[i] > 0
                cms[i] = cms[imax] #TODO: Probably not quite right, but shouldn't be an issue if stall not reached.
            else
                cms[i] = cms[imin]
            end
        end
    end
    

    return cls,cds,cms, clmax, clmin
end