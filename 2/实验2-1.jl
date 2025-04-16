#实验2，迭代法求解线性方程组,并判断收敛性
e=1e-6;#精度
Nmax=666; #迭代次数上限
A0=[2 1 0 4;3 7 -1 5;1 -1 4 4;1 2 3 -6];  #方程组(2-1)
b0=[10 12 15 3];
n=length(b0);

#Jacobi迭代法求解线性方程组，分量形式,(有2处需要根据注释填写代码)
x0=zeros(n,1);y=zeros(n,1);
x=x0;x0=x.+2*e;k=0;
b=b0';A=A0;
while norm(x0-x,Inf)>e&&k<Nmax
    global k,x,x0;
    k=k+1;x0=x;
    for i=1:n
        y[i]=b[i];
        for j=1:n
            if j != i
                y[i]=y[i]-A[i,j]*x0[j];      #(1)教材188页，公式（2.3）部分功能
            end
        end
        if abs(A[i,i])<1e-10
            warning("A (i, i) is too small, Jacobi iteration method fails");
            return
        end
       y[i]=y[i]/A[i,i];                #(2)实现教材188页，公式（2.3）部分功能
    end
    x=y;
end

if k==Nmax 
    println("Jacobi iteration method has failed,the maximum number of iterations has been reached.");
    println("  ");
else
   println("Jacobi Iterative method:");
   println("        solution vector              iterations\n");
   println("        $(round.(x',digits=3))                  $k");
end

#Gauss--Saidel迭代法求解线性方程组，矩阵形式(有1处需要根据注释填写代码)
x0=zeros(n,1);
x=x0;x0=x.+2*e;k=0;
b=b0';A=A0;
A1=tril(A);#提取矩阵下三角矩阵，对应D-L
iA1=inv(A1);#求逆矩阵
while norm(x0-x,Inf)>e&&k<Nmax
    global k,x,x0;
    k=k+1;x0=x;
    
    x=-iA1*(A-A1)*x0+iA1*b;  #(3)迭代计算解向量x，教材188页，公式（2.4）
   
end
if k==Nmax 
    println("Gauss--Saidel method has failed,the maximum number of iterations has been reached.");
else
    println("Gauss--Saidel Iterative method:");
    println("        solution vector              iterations\n");
    println("        $(round.(x',digits=3))                  $k");
end

# 两种迭代法的收敛性判定(有1处需要根据注释填写代码)
D=diagm(diag(A0));
J = I - inv(D) * A0 # 计算雅克比迭代矩阵
eigenvalues_J = eigen(J).values  # 提取特征值
pJ = maximum(abs.(eigenvalues_J))

println("Spectral radius of Jacobian iterative matrix:\n",pJ);
G=iA1*(triu(A0)-D);
eigenvalues_G = eigen(G).values  # 提取特征值
pG = maximum(abs.(eigenvalues_G))
println("Spectral radius of Gauss Saidel iterative matrix:\n",pG);