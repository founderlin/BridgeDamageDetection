# data analysis
ncores = 20
err_list = dataANY.err_cal(N, ncores)
# display(err_list)
err_df = convert(DataFrame, err_list)
# err_df = sort!(err_df, [:x4], rev = false)
err_df_s = sort!(err_df, [:x5], rev = false)
# err_df_s = sort!(err_df_s, [:x5], rev = false)
n_chosen = 100
err_df_sc = first(err_df_s, n_chosen)
err_ls_sc = convert(Matrix, err_df_sc)

i_err_ls = zeros(n_chosen, n+5)
for i = 1:n_chosen
   i_err = floor(Int, err_ls_sc[i,1])
   i_err_ls[i,1] = i_err_ls[i,1] + i_err
   i_err_ls[i,2:n+1] = i_err_ls[i,2:n+1] + ParaV[i_err,:]
   i_err_ls[i,n+2] = i_err_ls[i,n+2] + err_df_sc[i,2]
   i_err_ls[i,n+3] = i_err_ls[i,n+3] + err_df_sc[i,3]
   i_err_ls[i,n+4] = i_err_ls[i,n+4] + err_df_sc[i,4]
   i_err_ls[i,n+5] = i_err_ls[i,n+5] + err_df_sc[i,5]
end
i_err_df = convert(DataFrame, i_err_ls)
writedlm("results/rank_x5_1.txt", i_err_ls)
display(i_err_df)
