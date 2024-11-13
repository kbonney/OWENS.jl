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
import FLOWMath

path,_ = splitdir(@__FILE__)

# R_root = 0.001 # m biwing radius
# R_biwing = 20.0*3 # outer radius
# R_tip = 40.0*3 # outer radius
# nbelem_root = 5 #biwing elements for each 
# nbelem_biwing = 5 #tip elements
# nbelem_tip = 5 #tip elements

# N_control = 5
# bshapex_root = LinRange(0.0,R_root,N_control) #Blade shape, magnitude is relevant
# bshapez_root = zeros(N_control) #Blade shape, magnitude is relevant
# bshapex_biwing_U = LinRange(R_root,R_biwing,N_control) #Blade shape, magnitude is relevant
# bshapez_biwing_U = [0.0,2.0,2.0,1.0,0.0] #Blade shape, magnitude is relevant
# bshapex_biwing_L = LinRange(R_root,R_biwing,N_control) #Blade shape, magnitude is relevant
# bshapez_biwing_L = [-12.0,-8.0,-6.0,-4.0,0.0] #Blade shape, magnitude is relevant
# bshapex_tip = LinRange(R_biwing,R_tip,N_control) #Blade shape, magnitude is relevant
# bshapez_tip = zeros(N_control) #Blade shape, magnitude is relevant

# bladelen = sum(sqrt.((shapeX[2:end].-shapeX[1:end-1]).^2 .+ (shapeY[2:end].-shapeY[1:end-1]).^2 ))
# println("bladelen: $bladelen")


ntelem = 20 #tower elements
nbelem = 60 #blade elements
nselem = 10

delta_t=0.05
numTS=10
Nbld = 3
Blade_Height = 5.0
Blade_Radius = 100.0
RPM= 10.0
Vinf = 10.0
eta = 0.5
AModel="AD"
structuralModel = "GX"

shapeZ = [0,0.45,0.89,0.9,1].*Blade_Height
shapeX = [0.0,0.125,0.25,0.26,1].*Blade_Radius

spline_Z = LinRange(shapeZ[1],shapeZ[end],100)
spline_X = FLOWMath.akima(shapeZ,shapeX,spline_Z)

PyPlot.figure()
PyPlot.plot(shapeX,shapeZ,".-")
PyPlot.plot(spline_X,spline_Z,".-")

mymesh,myel,myort,myjoint,sectionPropsArray,mass_twr, mass_bld,
stiff_twr, stiff_bld,bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,aeroForces,deformAero,
mass_breakout_blds,mass_breakout_twr,system, assembly, sections,AD15bldNdIdxRng, AD15bldElIdxRng = OWENS.setupOWENS(OWENSAero,path;
    rho=1.225,
    Nslices=10,
    ntheta=30,
    RPM,
    Vinf,
    eta,
    B = Nbld,
    H = Blade_Height,
    R = Blade_Radius,
    shapeZ,
    shapeX,
    shapeY = zero(shapeX),
    ifw=true,
    delta_t,
    numTS,
    adi_lib=nothing,
    adi_rootname="$path/biwing",
    AD15hubR = 1.0,
    windINPfilename="$(path)/data/turbsim/350mx350m_30x30_20msETM.bts",
    ifw_libfile=nothing,
    NuMad_geom_xlscsv_file_twr = "$path/data/NuMAD_Geom_ARCUS330m_tower_DecDesign_noprebend_biwing.csv",
    NuMad_mat_xlscsv_file_twr = "$path/data/NuMAD_Materials_ARCUS330m_DecDesign.csv",
    NuMad_geom_xlscsv_file_bld = "$path/data/NuMAD_Geom_ARCUS330m_blades_DecDesign_noprebend_biwing.csv",
    NuMad_mat_xlscsv_file_bld = "$path/data/NuMAD_Materials_ARCUS330m_DecDesign.csv",
    NuMad_geom_xlscsv_file_strut = ["$path/data/NuMAD_Geom_SNL_5MW_Struts.csv"],
    NuMad_mat_xlscsv_file_strut = "$path/data/NuMAD_Materials_ARCUS330m_DecDesign.csv",
    Htwr_base=10.0,
    ntelem, 
    nbelem, 
    ncelem=5,
    nselem,
    joint_type = 0,
    strut_twr_mountpoint = [0.89],
    strut_bld_mountpoint = [0.89],
    AModel, #AD, DMS, AC
    DSModel="BV",
    RPI=true,
    cables_connected_to_blade_base = true,
    angularOffset = 0.0,
    meshtype = "Darrieus",
    isHAWT=true)


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
#### AeroElastic
#######################################
# node, dof, bc
pBC = [1 1 0
1 2 0
1 3 0
1 4 0
1 5 0
1 6 0]

if AModel=="AD"
    AD15On = true
else
    AD15On = false
end

tocp = [0.0,100000.1]
Omegaocp = [RPM,RPM] ./ 60
tocp_Vinf = [0.0,100000.1]
Vinfocp = [Vinf,Vinf]

inputs = OWENS.Inputs(;analysisType = structuralModel,
tocp,
Omegaocp,
tocp_Vinf,
Vinfocp,
numTS,
delta_t,
AD15On,
aeroLoadsOn = 2)

nothing

# Then there are inputs for the finite element models, also, please see the api reference for specifics on the options (TODO: ensure that this is propogated to the docs)

feamodel = OWENS.FEAModel(;analysisType = structuralModel,
outFilename = "none",
joint = myjoint,
platformTurbineConnectionNodeNumber = 1,
pBC,
nlOn = false,
gravityOn = [0,-9.81,0],
numNodes = mymesh.numNodes,
RayleighAlpha = 0.05,
RayleighBeta = 0.05,
iterationType = "DI")

nothing

# Here is where we actually call the unsteady simulation and where owens pulls the aero and structural solutions together
# and propogates things in time.

println("Running Unsteady")
t, aziHist,OmegaHist,OmegaDotHist,gbHist,gbDotHist,gbDotDotHist,FReactionHist,
FTwrBsHist,genTorque,genPower,torqueDriveShaft,uHist,uHist_prp,epsilon_x_hist,epsilon_y_hist,
epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,FPtfmHist,FHydroHist,FMooringHist,
topFexternal_hist,rbDataHist = OWENS.Unsteady_Land(inputs;system,assembly,
topModel=feamodel,topMesh=mymesh,topEl=myel,aero=aeroForces,deformAero)

azi=aziHist#./aziHist*1e-6
saveName = "$path/vtk/biwing"
tsave_idx=1:1:numTS-1
OWENS.OWENSVTK(saveName,t,uHist,system,assembly,sections,aziHist,mymesh,myel,
    epsilon_x_hist,epsilon_y_hist,epsilon_z_hist,kappa_x_hist,kappa_y_hist,kappa_z_hist,
    FReactionHist,topFexternal_hist;tsave_idx)

massOwens,stress_U,SF_ult_U,SF_buck_U,stress_L,SF_ult_L,SF_buck_L,stress_TU,SF_ult_TU,
SF_buck_TU,stress_TL,SF_ult_TL,SF_buck_TL,topstrainout_blade_U,topstrainout_blade_L,
topstrainout_tower_U,topstrainout_tower_L,topDamage_blade_U,
topDamage_blade_L,topDamage_tower_U,topDamage_tower_L = OWENS.extractSF(bld_precompinput,
bld_precompoutput,plyprops_bld,numadIn_bld,lam_U_bld,lam_L_bld,
twr_precompinput,twr_precompoutput,plyprops_twr,numadIn_twr,lam_U_twr,lam_L_twr,
mymesh,myel,myort,number_of_blades,epsilon_x_hist,kappa_y_hist,kappa_z_hist,epsilon_z_hist,
kappa_x_hist,epsilon_y_hist;verbosity, #Verbosity 0:no printing, 1: summary, 2: summary and spanwise worst safety factor # epsilon_x_hist_1,kappa_y_hist_1,kappa_z_hist_1,epsilon_z_hist_1,kappa_x_hist_1,epsilon_y_hist_1,
LE_U_idx=1,TE_U_idx=6,SparCapU_idx=3,ForePanelU_idx=2,AftPanelU_idx=5,
LE_L_idx=1,TE_L_idx=6,SparCapL_idx=3,ForePanelL_idx=2,AftPanelL_idx=5,
Twr_LE_U_idx=1,Twr_LE_L_idx=1,
AD15bldNdIdxRng,AD15bldElIdxRng,strut_precompoutput=nothing) #TODO: add in ability to have material safety factors and load safety factors


##############################################
# Modal
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
