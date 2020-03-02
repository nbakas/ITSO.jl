


import numpy as np
import copy

func_evals=10000;n_vars=10;
lb=-10*np.ones(n_vars);ub=15*np.ones(n_vars);
x=lb+(ub-lb)*np.random.random(n_vars);
x_evals=np.zeros((func_evals,n_vars));
x_evals[0,:]=lb
x_evals[1,:]=x
x_evals[2,:]=ub

def func(x):
    # return sum(x**2)
    return 1 + (1/4000)*sum((x+7)**2) - np.prod(np.cos((x+7) / range(1, len(x)+1, 1)))

opti_evals=np.zeros(func_evals)
l=len(opti_evals)
opti_evals[0]=func(lb)
opti_evals[1]=func(x)
opti_evals[2]=func(ub)
iopt=np.argsort(opti_evals[0:3])
opti_f=opti_evals[iopt[0]]
opti_x=x_evals[iopt[0],:]
x=copy.copy(opti_x)
iter=3

while iter <= func_evals:
    for j in range(0, n_vars, 1):
        iter+=1    
        if iter > func_evals:
            break
        if iter>func_evals/3: 
            alpha1=np.floor((func_evals/3)*(1-iter/func_evals))+2 
        else: 
            alpha1=iter 
        alpha1=int(alpha1)
        indsBEST=np.argsort(opti_evals[0:iter-1])[0:alpha1-1]
        inp_x=x_evals[indsBEST,j]
        rr=min(inp_x)+np.random.random(1)*(max(inp_x)-min(inp_x))
        x[j]=rr
        fu=func(x) 
        opti_evals[iter-1]=fu
        x_evals[iter-1,:]=copy.copy(x)

        if fu<opti_f:
            opti_x=copy.copy(x)
            opti_f=copy.copy(fu)
            print(" iter= ",iter," otpi_f= ",opti_f," alpha1= ",alpha1)
        else:
            x=copy.copy(opti_x)



print(" otpi_f= ",opti_f)
