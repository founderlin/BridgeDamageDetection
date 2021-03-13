@time ParaV, DF = paraVAR.para_var(n, rn)
# output = "para.csv"
# CSV.write(output, DF)

@time Um = moniSIM.static_cal(ma_para, f, n, ParaM)
display(Um)

# @time paraSIM.parallel_computation(N, n, rn, Um, ParaV, ma_para, f)

include("./module/paraSIM_X.jl")
