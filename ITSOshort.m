
clear all;
clc;

func_evals=10000;n_vars=10;
lb=-10*ones(n_vars,1);ub=15*ones(n_vars,1);
x=lb+(ub-lb).*rand(n_vars,1);
x_evals=zeros(func_evals,n_vars);
x_evals(1,:)=lb;x_evals(2,:)=x;x_evals(3,:)=ub;
function ret=func(x)
    # ret=sum(x.^2);
    # x=x';  # for ga
    n = length(x);ret=1 + (1/4000)*sum((x.+7).^2) - prod(cos((x.+7) ./ sqrt(1:n)'));
endfunction
opti_evals=zeros(func_evals,1);opti_evals(1)=func(lb);opti_evals(2)=func(x);opti_evals(3)=func(ub);
[so,iopt]=sort(opti_evals(1:3));opti_f=opti_evals(iopt(1));opti_x=x_evals(iopt(1),:)';
x=opti_x;iter=3;

while iter<func_evals
    for j=randi(n_vars)  %1:n_vars 
        iter+=1; if iter>func_evals return end
        if iter>func_evals/3 alpha1=floor((func_evals/3)*(1-iter/func_evals))+2; else alpha1=iter; end
        [so,indsBEST]=sort(opti_evals(1:iter-1));
        inp_x=x_evals(indsBEST(1:alpha1-1),j);
        rr=min(inp_x)+rand()*(max(inp_x)-min(inp_x));
        x(j)=rr; fu=func(x); opti_evals(iter)=fu; x_evals(iter,:)=x;

        if fu<opti_f opti_x=x;opti_f=fu;printf("iter= %d otpi_f= %f alpha1= %d \n", iter,opti_f,alpha1)
        else x=opti_x; end
    end
end
opti_evals.-=min(opti_evals); opti_evals./=max(opti_evals);
%plot(sort(opti_evals,'descend').^(1/10));
printf("iter= %d otpi_f= %f alpha1= %d \n", iter,opti_f,alpha1)


%options = gaoptimset('Generations', func_evals)
%[opti_x_ga, opti_f_ga, exit_flag]= ga(@func, n_vars, [], [], [], [], lb, ub, [], options)




