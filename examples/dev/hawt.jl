import PyPlot
PyPlot.pygui(true)
PyPlot.close("all")
PyPlot.rc("figure", figsize=(4.5, 3))
PyPlot.rc("font", size=10.0)
PyPlot.rc("lines", linewidth=1.5)
PyPlot.rc("lines", markersize=4.0)
PyPlot.rc("legend", frameon=true)
PyPlot.rc("axes.spines", right=false, top=false)
PyPlot.rc("figure.subplot", left=.18, bottom=.17, top=0.9, right=.9)
PyPlot.rc("figure",max_open_warning=500)
# PyPlot.rc("axes", prop_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
plot_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7"]

import OWENS
import OWENSAero
import GXBeam

path,_ = splitdir(@__FILE__)

# include("$(path)/../../../../OWENS.jl/src/OWENS.jl")
# include("$(path)/../src/OWENS.jl")
# include("$(path)/setupOWENShawt.jl")
println("Set up Macro Geometry/Inputs")
rho = 1.225
Nslices = 30
ntheta = 30
RPM = 7.8
Vinf = 20.0
eta = 0.5
B = Nbld = 3
R = 110.51 #m
H = 1.02*R*2 #m

ntelem = 20 #tower elements
nbelem = 60 #blade elements
ncelem = 10


shapeX = LinRange(0,R,Nslices+1)
shapeY = zero(shapeX)

R_root = 0.001 # m biwing radius
R_biwing = 20.0*3 # outer radius
R_tip = 40.0*3 # outer radius
nbelem_root = 5 #biwing elements for each 
nbelem_biwing = 5 #tip elements
nbelem_tip = 5 #tip elements

N_control = 5
bshapex_root = LinRange(0.0,R_root,N_control) #Blade shape, magnitude is relevant
bshapez_root = zeros(N_control) #Blade shape, magnitude is relevant
bshapex_biwing_U = LinRange(R_root,R_biwing,N_control) #Blade shape, magnitude is relevant
bshapez_biwing_U = [0.0,2.0,2.0,1.0,0.0] #Blade shape, magnitude is relevant
bshapex_biwing_L = LinRange(R_root,R_biwing,N_control) #Blade shape, magnitude is relevant
bshapez_biwing_L = [-12.0,-8.0,-6.0,-4.0,0.0] #Blade shape, magnitude is relevant
bshapex_tip = LinRange(R_biwing,R_tip,N_control) #Blade shape, magnitude is relevant
bshapez_tip = zeros(N_control) #Blade shape, magnitude is relevant

bladelen = sum(sqrt.((shapeX[2:end].-shapeX[1:end-1]).^2 .+ (shapeY[2:end].-shapeY[1:end-1]).^2 ))
println("bladelen: $bladelen")


mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
    stiff_twr, stiff_bld,bld_precompinput,
    bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
    twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,
    mass_breakout_blds,mass_breakout_twr,bladeIdx,bladeElem,system,assembly,sections = OWENS.setupOWENShawt(OWENSAero,path;
    rho,
    Nslices,
    ntheta,
    RPM,
    Vinf,
    eta,
    B = Nbld,
    H = 5.0,
    R = 2.5,
    hubR = 2.0,
    shapeZ = shapeY,
    shapeX,#shapeX_spline(shapeZ)
    ifw=false,
    windINPfilename="$(path)/data/turbsim/350mx350m_30x30_20msETM.bts",
    ifw_libfile = nothing,
    NuMad_geom_xlscsv_file_twr = "$path/data/NuMAD_Geom_ARCUS330m_tower_DecDesign_noprebend_biwing.csv",
    NuMad_mat_xlscsv_file_twr = "$path/data/NuMAD_Materials_ARCUS330m_DecDesign.csv",
    NuMad_geom_xlscsv_file_bld = "$path/data/NuMAD_Geom_ARCUS330m_blades_DecDesign_noprebend_biwing.csv",
    NuMad_mat_xlscsv_file_bld = "$path/data/NuMAD_Materials_ARCUS330m_DecDesign.csv",
    Htwr_base=10.0,
    ntelem, #tower elements
    nbelem, #blade elements
    ncelem,
    stack_layers_scale = [4.0,4.0],
    joint_type = 0,
    c_mount_ratio = 0.05,
    AModel="DMS",
    DSModel="BV",
    RPI=true,
    biwing=true,
    hub_depth = 15.0, #Hub Beam Depth
    R_root,
    R_biwing,
    R_tip,
    nbelem_root,
    nbelem_biwing,
    nbelem_tip,
    bshapex_root,
    bshapez_root,
    bshapex_biwing_U,
    bshapez_biwing_U,
    bshapex_biwing_L,
    bshapez_biwing_L,
    bshapex_tip,
    bshapez_tip,
    )


PyPlot.figure()
for icon = 1:length(mymesh.conn[:,1])
    idx1 = mymesh.conn[icon,1]
    idx2 = mymesh.conn[icon,2]
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"k.-")
    PyPlot.text3D(mymesh.x[idx1].+rand()/30,mymesh.y[idx1].+rand()/30,mymesh.z[idx1].+rand()/30,"$idx1",ha="center",va="center")
    # sleep(0.1)
end

for ijoint = 1:length(myjoint[:,1])
    idx2 = Int(myjoint[ijoint,2])
    idx1 = Int(myjoint[ijoint,3])
    PyPlot.plot3D([mymesh.x[idx1],mymesh.x[idx2]],[mymesh.y[idx1],mymesh.y[idx2]],[mymesh.z[idx1],mymesh.z[idx2]],"r.-")
    sleep(0.1)
end
    

println("\nBlades' Mass Breakout")
for (i,name) in enumerate(plyprops_bld.names)
    println("$name $(mass_breakout_blds[i]) kg, $(plyprops_bld.costs[i]) \$/kg: \$$(mass_breakout_blds[i]*plyprops_bld.costs[i])")
end

println("\nTower Mass Breakout")
for (i,name) in enumerate(plyprops_twr.names)
    println("$name $(mass_breakout_twr[i]) kg, $(plyprops_twr.costs[i]) \$/kg: \$$(mass_breakout_twr[i]*plyprops_twr.costs[i])")
end

println("Total Cost Blades: \$$(sum(mass_breakout_blds.*plyprops_bld.costs))")
println("Total Cost Tower: \$$(sum(mass_breakout_twr.*plyprops_twr.costs))")
println("Total Cost: \$$(sum(mass_breakout_blds.*plyprops_bld.costs)+ sum(mass_breakout_twr.*plyprops_twr.costs))")

println("\nBlades' Material Max Strain")
for (i,name) in enumerate(plyprops_bld.names)
    println("$name $(plyprops_bld.plies[i].xt/plyprops_bld.plies[i].e1) xt $(plyprops_bld.plies[i].xc/plyprops_bld.plies[i].e1) xc $(plyprops_bld.plies[i].yt/plyprops_bld.plies[i].e2) yt $(plyprops_bld.plies[i].yc/plyprops_bld.plies[i].e2) yc")
end

######################################
#### Perform Aerostructural One Way Test
#######################################

println("Saving VTK time domain files")
userPointNames=["EA","EIyy","EIzz","Fx","Fy","Fz","Mx","My","Mz"]
t = [0,1]

# map el props to points using con
userPointData = zeros(length(userPointNames),length(t),mymesh.numNodes)
EA_points = zeros(mymesh.numNodes)
EIyy_points = zeros(mymesh.numNodes)
EIzz_points = zeros(mymesh.numNodes)

# # Time-invariant data
# for iel = 1:length(myel.props)
#     nodes = mymesh.conn[iel,:]
#     EA_points[Int.(nodes)] = myel.props[iel].EA
#     EIyy_points[Int.(nodes)] = myel.props[iel].EIyy
#     EIzz_points[Int.(nodes)] = myel.props[iel].EIzz
# end

# # fill in the big matrix
# for it = 1:length(t)

#     userPointData[1,it,:] = EA_points
#     userPointData[2,it,:] = EIyy_points
#     userPointData[3,it,:] = EIzz_points
#     # userPointData[4,it,:] = FReactionHist[it,1:6:end]
#     # userPointData[5,it,:] = FReactionHist[it,2:6:end]
#     # userPointData[6,it,:] = FReactionHist[it,3:6:end]
#     # userPointData[7,it,:] = FReactionHist[it,4:6:end]
#     # userPointData[8,it,:] = FReactionHist[it,5:6:end]
#     # userPointData[9,it,:] = FReactionHist[it,6:6:end]
# end

azi=[0.0,pi/8]#./aziHist*1e-6
uHist = [zeros(mymesh.numNodes*6) zeros(mymesh.numNodes*6)]'
saveName = "$path/vtk/HAWTBiwingpleasework"
OWENS.OWENSFEA_VTK(saveName,t,uHist,system,assembly,sections;scaling=1,azi,userPointNames,userPointData)


# node, dof, bc
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

# filename = "$(path)/data/newmesh_34m"
# OWENS.saveOWENSfiles(filename,mymesh,myort,myjoint,myel,pBC,numadIn_bld)


# if testModal
##############################################
# Modal Test
#############################################
displ = zeros(mymesh.numNodes*6)
numModes = 32

FEAinputs = OWENS.FEAModel(;analysisType = "M",
        # outFilename = "$path/data/outplat.out",
        joint = myjoint,
        platformTurbineConnectionNodeNumber = 1,
        pBC = pBC,
        gravityOn=true,
        nlOn = true,
        tolerance = 1e-6,
        spinUpOn = true,
        iterationType = "NR",
        numNodes = mymesh.numNodes,
        numModes)  # number of modes to calculate)


starttime2 = time()
FEAinputs.analysisType = "GX"
freq2 = OWENS.AutoCampbellDiagram(FEAinputs,mymesh,myel,system,assembly,sections;
    minRPM = 0.0,
    maxRPM = 40.0,
    NRPM = 9, # int
    vtksavename="$path/campbellVTK/SNL34m",
    saveModes = [1,3,5], #must be int
    mode_scaling = 500.0,
    )
 
rotSpdArrayRPM = LinRange(0.0, 40.0, 9) # int 
freqGX = [freq2[:,i] for i=1:2:FEAinputs.numModes-6-2]
elapsedtime2 = time() - starttime2


PyPlot.figure()
# for i=1:1:numModes-2
#        PyPlot.plot(rotSpdArrayRPM,freqOW[:,i],color=plot_cycle[1],"b-") #plot mode i at various rotor speeds
# end

#plot per rev lines
for i=1:NperRevLines
    linex=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5]
    liney=[rotSpdArrayRPM[1], rotSpdArrayRPM[end]+5].*i./60.0
    PyPlot.plot(linex,liney,"--k",linewidth=0.5)
    PyPlot.annotate("$i P",xy=(0.95*linex[2],liney[2]+.05+(i-1)*.01))
end
PyPlot.grid()
PyPlot.xlabel("Rotor Speed (RPM)")
PyPlot.ylabel("Frequency (Hz)")
# PyPlot.plot(0,0,"k-",label="Experimental")
# PyPlot.plot(0,0,color=plot_cycle[1],"-",label="OWENS")
PyPlot.legend()
# PyPlot.ylim([0.0,0.8])
# PyPlot.savefig("$(path)/../figs/34mCampbell.pdf",transparent = true)

# Add to figure
for i=1:2:FEAinputs.numModes-6-2
       PyPlot.plot(rotSpdArrayRPM,freq2[:,i],color=plot_cycle[2],"-") #plot mode i at various rotor speeds
end
PyPlot.plot(0,0,color=plot_cycle[2],"-",label="GXBeam")
PyPlot.legend(fontsize=8.5,loc = (0.09,0.8),ncol=2,handleheight=1.8, labelspacing=0.03)
PyPlot.ylim([0,6.01])
# PyPlot.savefig("$(path)/../figs/34mCampbellWGX.pdf",transparent = true)
