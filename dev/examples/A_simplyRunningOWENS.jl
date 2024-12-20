import OWENS

runpath = path = "/home/runner/work/OWENS.jl/OWENS.jl/examples/literate" # to run locally, change to splitdir(@__FILE__)[1]

modelopt = OWENS.ModelingOptions("$(path)/OWENS_Opt.yml")
designparams = OWENS.Design_Data("$path/WINDIO_example.yaml")

OWENS.runOWENSWINDIO(modelopt,designparams,runpath)

modelopt.DLC_Options.DLCs = ["1_1"] #"normal"
#### modelopt.DLC_Options.DLCs = ["1_3","6_1"] #"normal"

OWENS.runDLC(modelopt,designparams,runpath)

nothing

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
