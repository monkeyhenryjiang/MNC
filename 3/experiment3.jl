#实验3，脚本文件，用不动点迭代、Newton法求解非线性方程
include("FixedPoint_Iteration.jl")
include("Newton_Iteration.jl")
# fun1=@(x)(x-2).^2;
# fplot(fun1,[0,4],'k','Linewidth',1.5);      #(x-2)^2
# fun2=@(x)log(x);
# fplot(fun2,[0,4],'r','Linewidth',1);        #lnx



x0=1:0.2:4;
n=length(x0);
result=zeros(2,n);              #初值对应的迭代次数统计
roots=zeros(2,n);               #初值对应的根
#不动点迭代法

fname=x->2+log(x).^0.5; #不动点迭代公式x=2+ln(x)^0.5

for i=1:n
    root,k=FixedPoint_Iteration(fname,x0[i],600,1e-8);
    result[1,i]=k;
    roots[1,i]=root;  
end

#牛顿迭代法
@syms x;        #符号变量，求导用
fun= x -> (x-2)^2-log(x); #原方程

for i=1:n
  if(x0[i] != 2.2)
     root,k=Newton_Iteration(fun,x0[i],600,1e-8);
     result[2,i]=k;
     roots[2,i]=root;
  end
     
end
 

#命令行窗口输出

fprintf("             Initial values:");   #迭代初值
 for i=1:n
     fprintf("%6.1f",x0[i]);
 end
 fprintf("\nFixed point iteration times:");  #不动点迭代次数
 for i=1:n
     fprintf("%6d",result[1,i]);
 end
 fprintf("\nFixed point iteration roots:");  
 for i=1:n
     fprintf("%6.3f",roots[1,i]);             #不动点迭代计算的根
 end
 fprintf("\n     Newton iteration times:");   #牛顿迭代次数
 for i=1:n
     fprintf("%6d",result[2,i]);
 end
 fprintf("\n     Newton iteration roots:");  
 for i=1:n
     fprintf("%6.3f",roots[2,i]);             #牛顿迭代计算的根
 end
 fprintf("\n");
#绘制直方图
data = hcat(result[1,:],result[2,:]);
ymax = maximum(result);
bar(x0,data,width = 1.2);
legend(["Fixed point iteration times","Newton iteration times"]);     #不动点迭代次数,牛顿迭代次数
xlabel("Initial values");ylabel("Iteration times");
title("The influence of initial value selection on the number of iterations");
