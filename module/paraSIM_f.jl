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
    dU = zeros(3*(n+1), sn)
    for i = 1:sn
        v = ParaV[i, :]
        K = fem.stiffMat(ea, ei, el, n, v);
        u = K^-1*f;
        du = u - um

        # Cross Entropy Loss
        # du_M = zeros(size(u))
        # for j = 1:size(u,2)
        #     for k = 1:size(u,1)
        #         println(um[k,j])
        #         du_item = -(u[k,j]*log(um[k,j]))+(1-u[k,j])*log(1-um[k,j]) # Cross Entropy Loss
        #         du_M[k,j] = du_item
        #     end
        # end

        # println(du)
        ddu = zeros(3*(n+1), 1)
        for j = 1:15
            norm_u = norm(du[j, :])
            ddu[j, 1] = norm_u # Square Error
            # ddu[j, 1] = sum(du[j, :])  # Bias Error
            # ddu[j, 1] = sum(du_M[j, :])  # Cross Entropy Loss-Bias Error
        end
        dU[:,i] = dU[:,i] + ddu
    end
    return dU
end

# function -> parallel simulation
function parallel_computation(N, n, rn, um, ParaV, ma_para, f, ncores::Int=21)
    dtn = floor(Int, N/(ncores-1))
    dUU = zeros(3*(n+1), N)
    @sync @distributed for i=1:ncores
        pa = 1 + dtn*(i-1)
        if dtn*i > N
            pb = N
            sn = N - pa + 1
        else
            pb = dtn*i
            sn = dtn
        end
        dU = parallel_static_cal(ma_para, f, sn, n, rn, um, ParaV[pa:pb, :])
        dUU[:, pa:pb] = dUU[:, pa:pb] + dU
        writedlm("./results/u_$(lpad(i,0)).txt", dUU'[:,[8,14]])
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
