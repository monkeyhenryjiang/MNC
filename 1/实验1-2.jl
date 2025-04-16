#实验1_2，分段三次(两点三次)Hermite插值法
hold("on");
fun(x)=5 / (2 + 5 * x^2);
x0=-6:0.01:6
y0=fun.(x0);
axis([-6.1 6.1 -0.5 3.5]);                      #设置坐标轴范围
plot(x0,y0,color=:black,linewidth=1);                                        #绘制原函数曲线
#分段三次Hermite插值法
x=-6:1:6;                               #插值点
n=length(x);
y=fun.(x);                                       
ydot=zeros(n);
x_i=-6:0.05:6;                              #绘图点
m=length(x_i);
y_i=zeros(m);                                   
@syms xx;                                       #符号变量，求导用
fname=5/(2 + xx^2*5);                             #符号方程，非函数
dfname=derivative(fname,xx);                              #derivative，求导数方程
df1=syslabFunction(dfname);                      #转换为函数句柄
ydot=df1.(x);
for i=1:m
    xi=x_i[i];
    yi=0;
    for k=1:n-1
      if x[k] <= xi && xi <= x[k+1]       #绘图点横坐标xi属于区间[x(k),x(k+1)]
      yi=y[k]*(1-2*(xi-x[k])/(x[k]-x[k+1]))*(xi-x[k+1])^2/(x[k]-x[k+1])^2+
      y[k+1]*(1-2*(xi-x[k+1])/(x[k+1]-x[k]))*(xi-x[k])^2/(x[k+1]-x[k])^2+
      ydot[k]*(xi-x[k])*(xi-x[k+1])^2/(x[k]-x[k+1])^2+
      +ydot[k+1]*(xi-x[k+1])*(xi-x[k])^2/(x[k+1]-x[k])^2;
      end
    end
    y_i[i]=yi;
end
ylabel("y=5/(2+5x^2)") ;          #设置坐标轴标注
plot(x_i,y_i,"r--",linewidth=2);
legend("实际曲线","Hermite插值曲线");
hold("off");
