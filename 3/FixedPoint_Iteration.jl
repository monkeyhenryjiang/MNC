#实验3_1，函数文件，不动点迭代求解非线性方程
function FixedPoint_Iteration(fname,x0,Nmax = 500,e = 1e-6)
# fname 为不动点迭代公式x=φ(x),x0为迭代初值；
# Nmax为最大的迭代次数，默认500；e为精度,默认值为1e-6
if fname === nothing || x0 === nothing
    error("Parameter error: fname and x0 are required.")
end

x=x0;x0=x+2*e;k=0;
while abs(x0-x)>e&&k<Nmax
    x0=x;
    x=fname(x0); #将x0带入到fun公式中，结果赋值给x
    k=k+1;
end
if x==Inf
    k=Nmax;root=Inf;
    return root,k;
end
if k==Nmax
    @warn("The calculation has failed，The maximum number of iterations has been reached.");
    root=Inf;
    return root,k;
else
    root=x;
    return root, k
end
end
