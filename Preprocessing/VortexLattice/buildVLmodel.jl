include("../../Structs/Aircraft.jl")
include("../Geometry.jl")
using VortexLattice

function surface_to_VL_surface(w::wing)
    xle = Vector{Float64}() # leading edge x-position
    yle = Vector{Float64}() # leading edge y-position
    zle = Vector{Float64}() # leading edge z-position
    chord = Vector{Float64}(w.chorddist) # chord length
    theta = Vector{Float64}(w.twistdist *pi/180) # twist (in radians)
    phi = zeros(Float64,length(w.points)) # section rotation about the x-axis
    fc = []#fill((xc) -> 0,length(w.points)) # camberline function for each section (y/c = f(x/c))
    for i = 1:1:length(w.points)
        push!(xle,w.points[i][1]+w.pos[1])
        push!(yle,w.points[i][2]+w.pos[2])
        push!(zle,w.points[i][3]+w.pos[3])
        sectionfc = airfoiltocamber(w.airfoils[i])
        push!(fc,sectionfc)#fc[i] = sectionfc
    end
    
    ns = 12 # number of spanwise panels
    nc = 6  # number of chordwise panels
    spacing_s = Uniform() # spanwise discretization scheme
    spacing_c = Uniform() # chordwise discretization scheme
    grid, surface = wing_to_surface_panels(xle, yle, zle, chord, theta, phi, ns, nc;fc = fc, spacing_s=spacing_s, spacing_c=spacing_c, mirror=false)
    global testxle = xle
    global testyle = yle
    global testzle = zle
    global testchord = chord
    global testtheta = theta
    global testphi = phi
    global testfc = fc
    global testsurf = surface
    global testref = Reference(1.0,1.0,1.0,zeros(3),1.0)
    global testfs = Freestream(1.0,0.001,0.0,[0.0;0.0;0.0])
    #global test = steady_analysis([surface],testref,testfs;symmetric = false)
    return surface
end

function buildVLmodel(craft::basicaircraft)
    VLsurfs = Vector{Matrix{SurfacePanel{Float64}}}()
    for surface in craft.surfs
        push!(VLsurfs,surface_to_VL_surface(surface))
    end
    #global testVLsurfs = VLsurfs #TEMP
    return VLsurfs
end