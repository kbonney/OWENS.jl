import OWENS

runpath = "/home/runner/work/OWENS.jl/OWENS.jl/docs/src/literate" #splitdir(@__FILE__)[1]

Inp = OWENS.MasterInput("/home/runner/work/OWENS.jl/OWENS.jl/docs/src/literate/sampleOWENS.yml")

OWENS.runOWENS(Inp,runpath)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
