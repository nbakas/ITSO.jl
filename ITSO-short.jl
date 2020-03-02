using Plots, Dates

func_evals=10000
n_vars=10
lb=-10ones(n_vars)
ub=15ones(n_vars)
x=lb+(ub-lb).*rand(n_vars)
x_evals=zeros(func_evals,n_vars)
x_evals[1,:]=lb
x_evals[2,:]=x
x_evals[3,:]=ub
function func(x)
    # sum(x.^2)

    n = length(x)
    1 + (1/4000)*sum(abs2, (x.+7)) - prod(cos.((x.+7) ./ sqrt.(1:n)))
end
opti_evals=zeros(func_evals)
opti_evals[1]=eval(func)(lb)
opti_evals[2]=eval(func)(x)
opti_evals[3]=eval(func)(ub)
iopt=sortperm(opti_evals[1:3])
opti_f=opti_evals[iopt[1]]
opti_x=x_evals[iopt[1],:]
x=copy(opti_x)
iter=3

while iter<func_evals
    global opti_f,opti_x,x_evals,opti_evals,iter
    for j=rand(1:n_vars)
        iter+=1; if iter>func_evals @goto escape_label end
        if iter>func_evals/3 alpha1=Int64(floor((func_evals/3)*(1-iter/func_evals)))+2 else alpha1=iter end
        # println(Dates.format(now(), "HH:MM:SS")," iter= ",iter," alpha1= ",alpha1)
        indsBEST=sortperm(opti_evals[1:iter-1])[1:alpha1-1]
        inp_x=x_evals[indsBEST,j];
        rr=minimum(inp_x)+rand()*(maximum(inp_x)-minimum(inp_x));
        x[j]=rr; fu=eval(func)(x); opti_evals[iter]=fu; x_evals[iter,:].=copy(x)

        if fu<opti_f
            opti_x.=copy(x);
            opti_f=copy(fu);
            println(Dates.format(now(), "HH:MM:SS")," iter= ",iter," otpi_f= ",opti_f," alpha1= ",alpha1)
        else
            x.=copy(vec(opti_x))
        end
    end
    @label escape_label
end
opti_evals.-=minimum(opti_evals); opti_evals./=maximum(opti_evals)
plot(sort(opti_evals,rev=true).^(1/10))



# using Optim
# res = Optim.optimize(func, fill(lb[1], n_vars), fill(ub[1], n_vars), x_evals[2,:], SAMIN(verbosity=0), Optim.Options(f_calls_limit=func_evals,store_trace =true));
# opti_x=res.minimizer; println(func(opti_x))
# res = Optim.optimize(func, fill(lb[1], n_vars), fill(ub[1], n_vars), x_evals[2,:], NelderMead(), Optim.Options(f_calls_limit=func_evals,store_trace =true))
# opti_x=res.minimizer; println(func(opti_x))
# res = Optim.optimize(func, fill(lb[1], n_vars), fill(ub[1], n_vars), x_evals[2,:], ParticleSwarm(), Optim.Options(f_calls_limit=func_evals,store_trace =true))
# opti_x=res.minimizer; println(func(opti_x))
