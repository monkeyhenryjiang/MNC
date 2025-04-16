#ʵ��3���ű��ļ����ò����������Newton���������Է���
include("FixedPoint_Iteration.jl")
include("Newton_Iteration.jl")
# fun1=@(x)(x-2).^2;
# fplot(fun1,[0,4],'k','Linewidth',1.5);      #(x-2)^2
# fun2=@(x)log(x);
# fplot(fun2,[0,4],'r','Linewidth',1);        #lnx



x0=1:0.2:4;
n=length(x0);
result=zeros(2,n);              #��ֵ��Ӧ�ĵ�������ͳ��
roots=zeros(2,n);               #��ֵ��Ӧ�ĸ�
#�����������

fname=x->2+log(x).^0.5; #�����������ʽx=2+ln(x)^0.5

for i=1:n
    root,k=FixedPoint_Iteration(fname,x0[i],600,1e-8);
    result[1,i]=k;
    roots[1,i]=root;  
end

#ţ�ٵ�����
@syms x;        #���ű���������
fun= x -> (x-2)^2-log(x); #ԭ����

for i=1:n
  if(x0[i] != 2.2)
     root,k=Newton_Iteration(fun,x0[i],600,1e-8);
     result[2,i]=k;
     roots[2,i]=root;
  end
     
end
 

#�����д������

fprintf("             Initial values:");   #������ֵ
 for i=1:n
     fprintf("%6.1f",x0[i]);
 end
 fprintf("\nFixed point iteration times:");  #�������������
 for i=1:n
     fprintf("%6d",result[1,i]);
 end
 fprintf("\nFixed point iteration roots:");  
 for i=1:n
     fprintf("%6.3f",roots[1,i]);             #�������������ĸ�
 end
 fprintf("\n     Newton iteration times:");   #ţ�ٵ�������
 for i=1:n
     fprintf("%6d",result[2,i]);
 end
 fprintf("\n     Newton iteration roots:");  
 for i=1:n
     fprintf("%6.3f",roots[2,i]);             #ţ�ٵ�������ĸ�
 end
 fprintf("\n");
#����ֱ��ͼ
data = hcat(result[1,:],result[2,:]);
ymax = maximum(result);
bar(x0,data,width = 1.2);
legend(["Fixed point iteration times","Newton iteration times"]);     #�������������,ţ�ٵ�������
xlabel("Initial values");ylabel("Iteration times");
title("The influence of initial value selection on the number of iterations");
