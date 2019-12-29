using Plots, BlackBoxOptim, Random, Dates, LinearAlgebra, Statistics, Printf, Optim, OptimTestProblems
plot(rand(10),xticks=1:1:10)

curr_path=pwd()
include(string(curr_path,"\\functions_opti.jl"))
include(string(curr_path,"\\call_bbox.jl"))
# ibb=[4 9 11] #
ibb=[1 2 4 5 6 7 9 10 11 12 13 14 15 16] #

funcs1=[:x_i, :elliptic, :cigar, :cigtab, :x_5_sq, :sin_01_x_2, :griewank, :quartic, :schwefel1_2, :rastrigin, :sphere, :ellipsoid, :alpine_1]
funcs1=[:sin_01_x_2]
# funcs1=[funcs1;funcs1;funcs1;funcs1;funcs1]

# include(string(curr_path,"\\dependencies_pdf_cdf.jl")) # για cdf
include(string(curr_path,"\\dependencies.jl")) # merged
vars1=Int64(20);evals1=Int64(10000)
history_all,times_all_step,x_evals,opti_evals,cdf_y,inp_yy,inp_x,pdf_y,iter,i_true=run_OHT(vars1,evals1,ibb,funcs1);
sizex=700; sizey=400; sf1=2.7
include(string(curr_path,"\\plots_all.jl"))


history_all_all=history_all.*0.0
nof_runs=10
for i=1:nof_runs
  global ibb,funcs1,history_all,vars1,evals1
  println("<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<",i,">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
  history_all,times_all_step,x_evals,opti_evals,cdf_y,inp_yy,inp_x,pdf_y,iter,i_true=run_OHT(vars1,evals1,ibb,funcs1);
  include(string(curr_path,"\\plots_all.jl"))
  history_all_all.+=history_all
end
history_all=history_all_all./nof_runs
sizex=500; sizey=400; sf1=2.7
include(string(curr_path,"\\plots_all.jl"))

savefig(string(curr_path,"\\FIGURES\\history\\all_",nof_runs,"times_",vars1,"vars_",evals1,"evals.pdf"))


