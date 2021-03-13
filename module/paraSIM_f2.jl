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

        # Normalisation
        u_ab = abs.(u)
        um_ab = abs.(um)
        u_n = zeros(size(u))
        um_n = zeros(size(um))
        for j = 1:size(u,2)
            u_n[:,j] = normalize(u_ab[:,j])
            um_n[:,j] = normalize(um_ab[:,j])
        end

        # Cross Entropy Loss
        du_M = zeros(size(u))
        for j = 1:size(u,2)
            for k = 1:size(u,1)
                du_item = -(u_n[k,j]*log(um_n[k,j]))+(1-u_n[k,j])*log(1-um_n[k,j]) # Cross Entropy Loss
                du_M[k,j] = du_item
            end
        end

        # println(du)
        ddu = zeros(3*(n+1), 1)
        for j = 1:15
            ddu[j, 1] = sum(du_M[j, :])  # Cross Entropy Loss-Bias Error
            # ddu[j, 1] = norm(du_M[j, :])  # Cross Entropy Loss-Square Error
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
