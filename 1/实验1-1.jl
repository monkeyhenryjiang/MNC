println("\033c");clear();                          #清屏和内存清除

fun(x)=5 ./ (2 .+ 5 .* x.^2);                       #定义函数f(x)
x_i=collect(range(-6, 6, length = 100));                       #取绘图点的横坐标值
p=length(x_i);                                   #计算绘图点数量

y_i_1=fun.(x_i);                             #计算原函数f(x)绘图的纵坐标值        

#牛顿插值法
x=-6:1.2:6;                                     #牛顿插值点，共11个
y=fun(x);  n=length(x);                                     
# 计算差商表
cssheet=zeros(n,n+1);
cssheet[:,1]=x';
cssheet[:,2]=y';                                   #差商表第2列

k=2;
while k !=n+1
    for i in k:n
     cssheet[i,k+1]=(cssheet[i,k]-cssheet[i-1,k])/(x[i]-x[i-k+1]);     #差商表计算
    end 
    global k = k + 1;                                   #此处要定义k是全局变量
end

# 计算newton插值
y_i_2=zeros(p);                               #定义牛顿插值函数绘图的纵坐标值   
for k = 1:p
    xi=x_i[k]; 
    cfwh=0;
    for i = 2:n
      w=1;
      for j = 1:i-1
          w=w*(xi-x[j]);
      end
    cfwh=cfwh+cssheet[i,i+1]*w;
    end
    y_i_2[k]=y[1]+cfwh;
end

#拉格朗日插值法
x=-6.1:0.6:6.1;                                #拉格朗日插值点，21个插值点
n=length(x);y=fun(x);
y_i_3=zeros(p);                                #定义拉格朗日插值绘图点纵坐标

for k=1:n
    t=ones(p);
    for j=1:n                                         #拉格朗日插值基函数
        if j != k
            t=t.*(x_i.-x[j])./(x[k]-x[j]);
        end
    end
    global y_i_3 = y_i_3 + t*y[k];
end  

#绘制3条曲线，依次是原函数曲线，牛顿插值曲线，拉格朗日插值曲线
hold("on");
plot(x_i,y_i_1,color=:black,linewidth=1.2);
plot(x_i,y_i_2,color=:blue,linewidth=1.5);
plot(x_i,y_i_3,"r-",linewidth=1.8);
xlim(-6.1,6.1),ylim(-4,10),title("插值方法");       #设置坐标抽范围，图题
ylabel("5/(5x^2+2)");
legend(["原函数" "牛顿插值曲线" "拉格朗日插值曲线"],loc="north"); hold("off");
