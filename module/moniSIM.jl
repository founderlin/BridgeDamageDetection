module moniSIM # monitoring data

include("./fem.jl")
# function -> static calculation
# --------------
    function static_cal(ma_para, f, n, ParaV)
        # ---------------------
        ea = ma_para[1];
        ei = ma_para[2];
        el = ma_para[3];

        # -----------------------
        v = ParaV
        K = fem.stiffMat(ea, ei, el, n, v);
        U = K^-1*f;
        return U
    end
end
