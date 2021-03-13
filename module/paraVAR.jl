module paraVAR # generate input parameter variants

using DataFrames

    function para_var(n, rn)

        # stochastic sampling
        # ----------------------------------------
        # v = rand(rn, n);
        v = [ 0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
              0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
              0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
              0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
              0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9;
              0.001 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
              ]'

        # Latin Hybercube Sampling
        # ----------------------------------------
        # v = lhs.LHS_all(rn, n)

        df = convert(DataFrame, v)
        df1 = df[:,[1]];
        for i = 2:n
             dfx = join(df1, df[:,[i]], kind = :cross)
             df1 = dfx
        end

        paraV = convert(Matrix, df1)
        ParaV = paraV
        # ParaV = zeros(rn^n, n*2)
        # for i = 1:n
        #     ParaV[:, 2*i-1] = df1[:, i]
        #     ParaV[:, 2*i] = df1[:, i]
        # end
        DF = convert(DataFrame, ParaV)
        # stochastic simulation
        N = size(ParaV , 1)

        # display(ParaV)
        return ParaV, DF
    end
end
