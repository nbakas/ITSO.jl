
BlackBoxOptim.MethodNames
BlackBoxNames=Array{String}(undef,0)
push!(BlackBoxNames,"Adapt. Diff. Evol. (rand 1 bin)")
push!(BlackBoxNames,"Adapt. Diff. Evol. (rand 1 bin, radious limited)")
push!(BlackBoxNames,"Borg Multiobj. Evol. Alg.")
push!(BlackBoxNames,"Diff. Evol. (rand 1 bin)")
push!(BlackBoxNames,"Diff. Evol. (rand 1 bin, radious limited)")
push!(BlackBoxNames,"Diff. Evol. (rand 2 bin)")
push!(BlackBoxNames,"Diff. Evol. (rand 2 bin, radious limited)")
push!(BlackBoxNames,"Distance-weighted Exponential Nat. Evol. Strat.")
push!(BlackBoxNames,"Compass Coordinate Direct Search")
push!(BlackBoxNames,"Probabilistic Descent Direct search")
push!(BlackBoxNames,"Random Search")
push!(BlackBoxNames,"Resampling Inheritance Memetic Search")
push!(BlackBoxNames,"Resampling Memetic Search")
push!(BlackBoxNames,"Separable Natural Evolution Strategies")
push!(BlackBoxNames,"Simultaneous Perturbation Stochastic Approximation")
push!(BlackBoxNames,"Exponential Natural Evolution Strategies")

ibb2=[ibb 17 18 19 20]
history_plot=history_all[1:iter,ibb2[1,:],:] # func_evals,n_optimizers,length(funcs)
n_optimizers=length(ibb2)
history_plotINIT=copy(history_plot)
for i=1:n_optimizers
    for j=1:size(history_plot,3)
        history_plotINIT[1:3,i,j].=maximum(history_plotINIT[1:3,:,j])
        mi=minimum(history_plotINIT[:,:,j])
        ma=maximum(history_plotINIT[:,:,j])
        # println(mi," ",ma)
        history_plot[:,i,j]=(history_plotINIT[:,i,j].-mi)./(ma.-mi)
    end
end
history_plot=mean(history_plot, dims=3)[:,:,1]
mi=minimum(history_plot[end,:])
ma=maximum(history_plot[end,:])
sR=history_plot[end,length(ibb)+1]
ssR=@sprintf("mi= %.2E ma= %.2E ITS= %.2E",mi,ma,sR)
println(ssR)
# display(plot(history_plot,width=1,label=[string.(ibb) "oht"],legend=:topleft,title=ssR))
display(plot(history_plot[:,1:length(ibb)],width=1,label=BlackBoxNames[ibb],legend=:topright,title=ssR,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+2],width=1,label="Nelder-Mead",legend=:topright,title=ssR,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+3],width=1,label="Particle Swarm",legend=:topright,title=ssR,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+4],width=1,label="Simulated Annealing",legend=:topright,title=ssR,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+1],width=2,label="Proposed-BNO",legendfontsize=8,color=:black,size = (sizex, sizey)))
display(xlims!(0,sf1*evals1))
display(xticks!(0:1000:evals1))
xaxis!("Function Evaluations")
yaxis!("Normalized Responces")
for i=[1:length(ibb);length(ibb)+2;length(ibb)+3;length(ibb)+4]
    println(minimum(history_plot[:,i]))
end
println(minimum(history_plot[:,length(ibb)+1]))

history_plot=history_all[1:iter,ibb2[1,:],:] # func_evals,n_optimizers,length(funcs)
n_optimizers=length(ibb2)
history_plotINIT=copy(history_plot)
for i=1:n_optimizers
    for j=1:size(history_plot,3)
        history_plotINIT[1:3,i,j].=maximum(history_plotINIT[1:3,:,j])
        mi=minimum(history_plotINIT[:,:,j])
        ma=maximum(history_plotINIT[:,:,j])
        history_plot[:,i,j]=(history_plotINIT[:,i,j].-mi)./(ma.-mi)
    end
end
# history_plot=(mean(history_plot, dims=3)[:,:,1].^(1/10))
history_plot=mean(history_plot, dims=3)[:,:,1].^(1/10)
# mi=minimum(history_plot[end,:])
# ma=maximum(history_plot[end,:])
# sR=history_plot[end,length(ibb)+1]
# ssR=@sprintf("mi= %.2E ma= %.2E ITS= %.2E",mi,ma,sR)
# println(ssR)
# display(plot(history_plot,width=1,label=[string.(ibb) "oht"],legend=:topleft,title=ssR))
display(plot(history_plot[:,1:length(ibb)],width=1,label=BlackBoxNames[ibb],legend=:false,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+2],width=1,label="Nelder-Mead",legend=:false,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+3],width=1,label="Particle Swarm",legend=:false,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+4],width=1,label="Simulated Annealing",legend=:false,titlefont = font("Times", 9),size = (sizex, sizey)))
display(plot!(history_plot[:,length(ibb)+1],width=2,label="Proposed-BNO",legend=:false,color=:black,size = (sizex, sizey)))
# display(xlims!(0,sf1*evals1))
display(xticks!(0:1000:evals1))
xaxis!("Function Evaluations")
yaxis!("Normalized Responces")
# xlims()


