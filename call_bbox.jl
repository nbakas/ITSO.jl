function call_bbox(bb_methods,all_bounds,func_evals,j,func,max_time)

    opt1 = bbsetup(@eval $func; Method=bb_methods[j], SearchRange = all_bounds, TraceMode = :silent, MaxFuncEvals = func_evals)# , MaxTime = max_time
    res1 = bboptimize(opt1)
    opti1=res1.archive_output.best_fitness
    sol1=res1.archive_output.best_candidate
    hist_num = length(opt1.runcontrollers[1].evaluator.archive.fitness_history)
    hist_tmp=zeros(hist_num)
    for i=1:hist_num
            hist_tmp[i]=opt1.runcontrollers[1].evaluator.archive.fitness_history[i].fitness
    end
    matt=repeat(hist_tmp', trunc(Int, func_evals/hist_num)+1, 1)[:]
    matt[end]=copy(opti1)

    return matt,opti1,sol1
end

#
