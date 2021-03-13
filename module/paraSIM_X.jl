# module paraSIM # execute parallel simulation

using Distributed

# function -> static calculation
# --------------
@everywhere using LinearAlgebra
@everywhere using DelimitedFiles
@everywhere include("./module/fem.jl")

@everywhere function parallel_static_cal(ma_para, f, sn, n, rn, um, ParaV)
    # ---------------------
    ea = ma_para[1];
    ei = ma_para[2];
    el = ma_para[3];

    # -----------------------
    dU = zeros(2, sn)
    dCt = zeros(1, sn)
    for i = 1:sn
        v = ParaV[i, :]
        K = fem.stiffMat(ea, ei, el, n, v);
        u = K^-1*f;

        cnt = 0
        u_c = u[[8,14],:]
        um_c = um[[8,14],:]

        du = (u_c - um_c)/um_c
        # println(du)
        ddu = zeros(2, 1)
        for j = 1:2
            norm_u = norm(du[j, :])
            ddu[j, 1] = norm_u # Square Error
            if norm_u < 0.1
                cnt = cnt+1
            end
        end
        dCt[i] = dCt[i] + cnt
        dU[:,i] = dU[:,i] + ddu
    end
    dU, dCt
end

# function -> parallel simulation
function parallel_computation(N, n, rn, um, ParaV, ma_para, f, ncores::Int=21)
    dtn = floor(Int, N/(ncores-1))
    dUU = zeros(2, N)
    dCT = zeros(1, N)
    @sync @distributed for i=1:ncores
        pa = 1 + dtn*(i-1)
        if dtn*i > N
            pb = N
            sn = N - pa + 1
        else
            pb = dtn*i
            sn = dtn
        end
        dU,dCt = parallel_static_cal(ma_para, f, sn, n, rn, um, ParaV[pa:pb, :])
        # dCt = parallel_match_cnt(ma_para, f, sn, n, rn, um, ParaV[pa:pb, :])
        dUU[:, pa:pb] = dUU[:, pa:pb] + dU
        dCT[:, pa:pb] = dCT[:, pa:pb] + dCt
        writedlm("./results/u_$(lpad(i,0)).txt", dUU')
        writedlm("./results/m_$(lpad(i,0)).txt", dCT')
    end
end
# end

# @time U = static_cal(N, n, rn, ParaV)
# display(U')
# U_df = convert(DataFrame, U')
# display(U_df)
# output = "u.csv"
# CSV.write(output, U_df)
@time parallel_computation(N, n, rn, Um, ParaV, ma_para, f)
# UU_df = convert(DataFrame, UU)
# display(UU_df)
# output = "uu.csv"
# CSV.write(output, UU_df)
