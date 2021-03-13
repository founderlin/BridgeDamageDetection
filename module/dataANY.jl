module dataANY

using DelimitedFiles
using LinearAlgebra

function collect(N, ncores)
    du = zeros(N, 2)
    dc = zeros(N, 1)
    for i=1:ncores
        ddu = readdlm("results/u_$(lpad(i,0)).txt")
        ddc = readdlm("results/m_$(lpad(i,0)).txt")
        du = du + ddu
        dc = ddc
    end
    return du, dc
end

function err_cal(N, ncores)
    du, dc = collect(N, ncores)
    err_list = zeros(N, 5)
    for i=1:N
        # err = norm(du[i,:])
        err = du[i,:].*10^3
        cnt = dc[i]
        err_list[i,1] = err_list[i,1] + i
        err_list[i,2:3] = err_list[2:3] + err
        err_list[i,4] = err_list[i,4] + abs(err[1]-err[2])
        # err_list[i,5] = err_list[i,5] + abs(err[1]+err[2])
        err_list[i,5] = err_list[i,5] + norm(err)
    end
    return err_list
end

end
