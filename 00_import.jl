# import packages
# --------------
# using CSV
using DelimitedFiles
using DataFrames
using Distributed

# import modules
# --------------
include("./module/fem.jl")
include("./module/paraVAR.jl")
include("./module/moniSIM.jl")
# include("./module/paraSIM_f.jl")
include("./module/dataANY.jl")


# add parallel proess
# --------------
# addprocs(21)
