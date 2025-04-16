using ForwardDiff
#实验3_2，函数文件，牛顿法求解非线性方程
function Newton_Iteration(fun,x0,Nmax = 500,e = 1e-6)
# 函数功能：返回根和迭代次数
#fun为原方程；Nmax为最大的迭代次数，默认500；e为精度,默认值为1e-6
if fun === nothing || x0 === nothing
    error("参数错误: fun 和 x0 是必需的。")
end
dfun(x) = ForwardDiff.derivative(fun, x);  #原方程求导
fname(x) = x - fun(x) / dfun(x); #牛顿迭代公式,隐函数形式
# fprintf('x0=#f\n',x0);
# fprintf('fun(x)=#f\n',fname(x0));
# fprintf('fun(x)=#f\n',dfun(x0));
xk=x0;  xkplus1=x0+2*e;   #初始化x(k)和x(k+1)
k=0;
while abs(xkplus1-xk) > e && k < Nmax
   if abs(dfun(xk))<e 
     fprintf("初值为%.1f时，出现导数值过小，",x0);
     @warn("停止计算");
     k = Nmax;  #出现导数值过小的特殊情况时直接将迭代次数置为500
     break;
   end
   if(k != 0)
    xk =  xkplus1;
   end
   xkplus1 = fname(xk);   #代入x(k)，计算x(k+1)
   k = k+1;
if k == Nmax 
    @warn("已达到迭代上限次数,停止计算"); 
      break;
 end
end

if k<Nmax && xkplus1==Inf
    k=Nmax;root=Inf;
    return root,k;
end
if k==Nmax
    @warn("已达到迭代次数上限，计算失败");
    root=Inf;
    return root,k;
else
    root = xkplus1;
    return root,k;
end
end