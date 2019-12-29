

function init_all_OHT(n_vars::Int64,func_evals::Int64,nof_funcs::Int64,nof_optimizers::Int64)
    println(Dates.format(now(), "HH:MM:SS")," >>>>>>>>>>>>>START----------------------------------------------->")
    bb_methods = copy(BlackBoxOptim.MethodNames)
    history_all=zeros(Float64,func_evals,17+3,nof_funcs).+1e10
    domain_width=10
    lb=-domain_width*ones(Float64,n_vars)
    ub=-2lb
    all_bounds=Tuple{Float64,Float64}[(lb[z], ub[z]) for z in 1:n_vars]
    optis=zeros(Float64,nof_optimizers+1+3)
    times_all=zeros(nof_optimizers+1+3)
    x_evals=zeros(func_evals,n_vars)
    opti_evals=zeros(func_evals).+1e100
    return history_all, lb, ub, all_bounds, optis, times_all, bb_methods, x_evals, opti_evals
end
function init_vars_OHT(func::Symbol,lb::Array{Float64,1},ub::Array{Float64,1},history_all::Array{Float64,3},n_vars::Int64,func_evals::Int64,fff::Int64,x_evals::Array{Float64,2},opti_evals::Array{Float64,1})
    x=lb+(ub-lb).*rand(n_vars)
    x_evals[1,:]=lb
    x_evals[2,:]=x
    x_evals[3,:]=ub
    opti_evals[1]=eval(func)(lb)
    opti_evals[2]=eval(func)(x)
    opti_evals[3]=eval(func)(ub)
    iopt=sortperm(opti_evals[1:3])
    opti_f=opti_evals[iopt[1]]
    opti_x=x_evals[iopt[1],:]
    x=copy(opti_x)
    iter=3
    history_all[1:3,17,fff].=opti_evals[1:3]
    return x, func, x_evals, opti_evals, opti_f, opti_x, iter, history_all
end
function run_OHT(n_vars::Int64,func_evals::Int64,ibb::Array{Int64,2},funcs::Array{Symbol,1})
    history_all, lb, ub, all_bounds, optis, times_all, bb_methods, x_evals, opti_evals = init_all_OHT(n_vars,func_evals,length(funcs),length(ibb))
    times_all_step=Array{Float64,1}(undef,0);iter=0;min_iter=Inf
    SR_all=zeros(length(funcs))
    alpha_plus_x=0; alpha_plus_f=0; alpha_minus_f=0; i_true=trues(func_evals)
    cdf_y=[];inp_yy=[];inp_y=[];inp_x=[];pdf_y=[];opti_f=Inf;indsBEST=[];var_xx=[];alpha1=[];var_yy=[];
    delta_fu_all=zeros(Float64,func_evals,n_vars)
    # delta_fu_all[1:3,:].=1.0
    for fff=1:length(funcs)
        nof_a1=[]
        alpha_plus_x=0; alpha_plus_f=0; alpha_minus_f=0; i_true=trues(func_evals)
        cdf_y=[];inp_yy=[];inp_y=[];inp_x=[];pdf_y=[];opti_f=Inf;indsBEST=[];var_xx=[];alpha1=[];var_yy=[];
        # try
            t1=now()
            x, func, x_evals, opti_evals, opti_f, opti_x, iter, history_all = init_vars_OHT(funcs[fff],lb,ub,history_all,n_vars,func_evals,fff,x_evals,opti_evals)
            alpha_plus_x=0; alpha_plus_f=0; alpha_minus_f=0; i_true=trues(func_evals,n_vars)
            alpha1s=[]
            while iter<func_evals
                for i=1:n_vars
                    for fght=1:1
                        alpha1=copy(iter)
                        if iter>func_evals/5 alpha1=Int64(floor((func_evals/5+10)*(1-iter/func_evals)^1)) +12 end
                        if alpha1>iter alpha1=copy(iter) end; push!(alpha1s,alpha1)
                        indsBEST=sortperm(opti_evals[1:iter])[1:alpha1]

                        if alpha1>50
                            curr_dists=x_evals[indsBEST,:].-repeat(opti_x',alpha1,1)
                            i_others=deleteat!(collect(1:n_vars),i)
                            curr_dists=sum(abs.(curr_dists[:,i_others]),dims=2)[:,1]
                            i_near=sortperm(curr_dists)[1:Int64(floor(alpha1/10))]
                            indsBEST=indsBEST[i_near]
                        end

                        inp_x=x_evals[indsBEST,i]; inp_y=opti_evals[indsBEST]; inp_d=delta_fu_all[indsBEST,i]
                        ind_sort_x=sortperm(inp_x); inp_x=inp_x[ind_sort_x]; inp_y=inp_y[ind_sort_x];

                        # if iter in [40;100;500;900] display(scatter(inp_x,inp_y,size=(320,300),legend=false,title=string(iter," ",i))) end

                        rr=rand()
                        if rr>1-iter/func_evals+0.1  # 0.5
                            x[i]=inp_x[1] .+ rand().*(inp_x[end].-inp_x[1])
                        else
                            x[i]=lb[i] .+ rand().*(ub[i].-lb[i])
                        end

                        iter+=1; x_evals[iter,:].=copy(x); fu=eval(func)(x);opti_evals[iter]=fu


                        if fu<opti_f opti_x.=copy(x); opti_f=copy(fu);  else x.=copy(vec(opti_x)) end
                        history_all[iter,17,fff]=opti_f;if iter==func_evals @goto  escape_label   end

                    end
                end
            end
            @label escape_label
            push!(times_all_step,convert(Int64,Dates.value(now()-t1))/1000)
            # println("unique f=",length(unique(opti_evals[1:iter]))," iter=",iter," opti_f=",opti_f," opti_x=",opti_x)
            if iter<func_evals history_all[iter:func_evals,17,fff].=opti_f end
            if iter<min_iter min_iter=copy(iter) end
            optis[end-3]=copy(opti_f)
            times_all[end-3]=convert(Int64,Dates.value(now()-t1))/1000
            # display(scatter(alpha1s,size=(320,300)))
            indt=Int64(0)
            for j=ibb
                indt+=1
                t1=now()
                matt,opti1,sol1 = call_bbox(bb_methods,all_bounds,iter,j,func,times_all[end])
                # println(sol1)
                history_all[1:iter,j,fff]=matt[1:iter]'
                optis[indt]=opti1
                times_all[indt]=convert(Int64,Dates.value(now()-t1))/1000
            end
            t1=now()
            res = Optim.optimize(getfield(Main, func), fill(lb[1], n_vars), fill(ub[1], n_vars), x_evals[2,:], NelderMead(), Optim.Options(f_calls_limit=func_evals,store_trace =true))
            # println(res.minimizer)
            hh=Optim.f_trace(res);hh[findfirst(hh.==minimum(hh)):end].=minimum(hh);matt=repeat(hh', trunc(Int, func_evals/length(hh))+1, 1)[:]
            history_all[1:iter,18,fff]=matt[1:iter]'
            optis[end-2]=res.minimum
            times_all[end-2]=convert(Int64,Dates.value(now()-t1))/1000

            res = Optim.optimize(getfield(Main, func), fill(lb[1], n_vars), fill(ub[1], n_vars), x_evals[2,:], ParticleSwarm(), Optim.Options(f_calls_limit=func_evals,store_trace =true))
            # println(res.minimizer)
            hh=Optim.f_trace(res);hh[findfirst(hh.==minimum(hh)):end].=minimum(hh);matt=repeat(hh', trunc(Int, func_evals/length(hh))+1, 1)[:];
            history_all[1:iter,19,fff]=matt[1:iter]'
            optis[end-1]=res.minimum
            times_all[end-1]=convert(Int64,Dates.value(now()-t1))/1000

            res = Optim.optimize(getfield(Main, func), fill(lb[1], n_vars), fill(ub[1], n_vars), x_evals[2,:], SAMIN(verbosity=0), Optim.Options(f_calls_limit=func_evals,store_trace =true));
            # println(res.minimizer)
			matt=[]
            try hh=Optim.f_trace(res);hh[findfirst(hh.==minimum(hh)):end].=minimum(hh);matt=repeat(hh', trunc(Int, func_evals/length(hh))+1, 1)[:]; catch
                matt=collect(range(opti_evals[2],length=func_evals,stop=res.minimum)); end;
            history_all[1:iter,20,fff]=matt[1:iter]';
            optis[end]=res.minimum;
            times_all[end]=convert(Int64,Dates.value(now()-t1))/1000


            SR_all[fff]=(opti_f-minimum(optis))/maximum(opti_evals)
            println(fff," of ",length(funcs)," ",Dates.format(now(), "HH:MM:SS")," func=", string(funcs[fff]), " ",optis," first_fu= ",convert(Int64,floor(opti_evals[1]))," times= ",times_all," SR_all[fff]= ",SR_all[fff],"\n")

            # catch ex
        #     if isa(ex, InterruptException) @warn("InterruptException..."); break
        #     else println(iter," ",ex," ",inp_y," ",pdf_y)
        #         if iter<min_iter min_iter=copy(iter) end end
        # end
    end
    # display(bar(SR_all,size=(320,300)))
    # println(delta_fu_all)
    return history_all, times_all_step, x_evals, opti_evals, cdf_y,inp_yy,inp_x,pdf_y,min_iter,i_true
end



# inp_yy[inp_yy.<1e-5].=0.0
# pdf_y=exp.(-(20-16iter/func_evals)inp_yy.^2)
# pdf_y=exp.(-4inp_yy.^2)
# pdf_y=1 ./(.1 .+inp_yy.^2)
# pdf_y=1 ./(inp_yy .+1e-100 .+(1e-10)*(1-iter/func_evals))
# pdf_y=1 ./(inp_yy .+1e-2)
# pdf_y=1 ./(inp_yy .+1e-2*(1-iter/func_evals))
