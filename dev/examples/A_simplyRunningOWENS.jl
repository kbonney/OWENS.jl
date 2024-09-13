import OWENS

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]
# runpath = path = splitdir(@__FILE__)[1]

Inp = OWENS.MasterInput("$runpath/sampleOWENS.yml")

OWENS.runOWENS(Inp,runpath)

simulated_time = 2.0 #seconds
DLCs = ["1_1"] #"normal"
#### DLCs = ["1_3"] #"normal"
#### DLCs = ["1_4"] #"normal"
#### DLCs = ["1_5"] #"normal"
#### DLCs = ["2_1"] #"freewheelatNormalOperatingRPM"
#### DLCs = ["2_3"] #"freewheelatNormalOperatingRPM"
#### DLCs = ["3_1"] #"startup"
#### DLCs = ["3_2"] #"startup"
#### DLCs = ["3_3"] #"startup"
#### DLCs = ["4_1"] #"shutdown"
#### DLCs = ["4_2"] #"shutdown"
#### DLCs = ["5_1"] #"emergencyshutdown"
#### DLCs = ["6_1"] #"parked"
#### DLCs = ["6_2"] #"parked_idle"
#### DLCs = ["6_4"] #"parked"
#### DLCs = ["7_1"] #"parked"
#### DLCs = ["2_3","3_1","3_2","3_3","4_1","4_2","5_1"]

OWENS.runDLC(DLCs,Inp,runpath;
    IEC_std="\"1-ED3\"",
    WindChar="\"A\"",
    WindClass=1,
    NumGrid_Z=38,
    NumGrid_Y=26,
    Vdesign=11.0,
    grid_oversize=1.25,
    Vinf_range=[10.0],#LinRange(4,24,21),
    regenWindFiles=true,
    delta_t_turbsim=0.05,
    simtime_turbsim=30.0,
    pathtoturbsim="$(OWENS.OWENSOpenFASTWrappers.OFWpath)/../deps/openfast/build/modules/turbsim/turbsim",
    runScript=OWENS.runOWENS)

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
