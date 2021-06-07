%LBM-MRT 2-D2Q9, note that c2=1/3, w9=4/9,
% 形式：f=zeros(9,nx,ny,);
% w1-4=1/9, and w5-w8, 1/36
clear
nx=101;ny=31;
f=zeros(9,nx,ny);feq=zeros(9,nx,ny); Fi=zeros(9,nx,ny);
u=zeros(nx,ny);v=zeros(nx,ny);
rho=ones(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
c2=1./3.;
L=ny-1;
uo=0.10;
Re=10;
alpha=uo*L/Re;
omega=1./(3.*alpha+0.5);
forcex=8.*uo*alpha./(L.^2);
forcey=0;
count=0; tol=1.0e-8; error=10.;erso=0.0;
s1=[omega omega-0.1 0 1.2 0 1.2 omega omega 0];
M=[1 1 1 1 1 1 1 1 1;-4 -1 -1 -1 -1 2 2 2 2;4 -2 -2 -2 -2 1 1 1 1;...
    0 1 0 -1 0 1 -1 -1 1;0 -2 0 2 0 1 -1 -1 1;0 0 1 0 -1 1 1 -1 -1;...
    0 0 -2 0 2 1 1 -1 -1;0 1 -1 1 -1 0 0 0 0;0 0 0 0 0 1 -1 1 -1;];
mm=circshift(M,[-1 -1]);
mminv=inv(mm);
msl=mminv*diag(s1);
%初始化
for i=1:nx
     for j=1:ny
          for k=1:9
             f(k,i,j)=rho(i,j)*w(k);%  计算初始态下的  分布函数=平衡分布函数                        
             Fi(k,i,j)=0;   %初始化离散的边界作用力分布函数,初始值赋为 0
         end
      end
end
 tic;
%Main Loop
while error>tol
    % Collitions
    [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,s1,forcex,forcey);
    % Streaming:
    [f]=stream(f,cx,cy);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,uo,rho);
    % Calculate rho, u, v
    [rho,u,v]=ruv(f,nx,ny,cx,cy,forcex,forcey);
    count=count+1;
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j);
        end
    end
    error=abs(ers-erso);
    erso=ers;
    if mod(count,50)==0
        result(nx,ny,u,v,count);
    end
end
result(nx,ny,u,v,count);
toc;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [f]=boundary(nx,ny,f,uo,rho)
%Boundary condition:
for j=1:ny
    f(3,nx,j)=f(3,2,j);
    f(7,nx,j)=f(7,2,j);
    f(6,nx,j)=f(6,2,j);
    f(4,1,j)=f(4,nx-1,j);
    f(7,1,j)=f(7,nx-1,j);
    f(8,1,j)=f(8,nx-1,j);  
end
% %bottom, and top boundary, bounce back
for i=1:nx
    f(2,i,1)=f(4,i,1);
    f(5,i,1)=f(7,i,1);
    f(6,i,1)=f(8,i,1);
    f(4,i,ny)=f(2,i,ny);
    f(7,i,ny)=f(5,i,ny);
    f(8,i,ny)=f(6,i,ny);
end

%Left boundary, velocity is given= uo
% for j=2:ny-1
%     f(3,nx,j)=f(3,nx-1,j);
%     f(6,nx,j)=f(6,nx-1,j);
%     f(7,nx,j)=f(7,nx-1,j);
%     f(1,1,j)=f(3,1,j)+2.*rho(1,j)*uo/3.;
%     f(5,1,j)=f(7,1,j)-0.5*(f(2,1,j)-f(4,1,j))+rho(1,j)*uo/6.;
%     f(8,1,j)=f(6,1,j)+0.5*(f(2,1,j)-f(4,1,j))+rho(1,j)*uo/6.;
% end
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Collition
function [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,s1,forcex,forcey)
% Calculate the
for i=1:nx
    for j=1:ny
        m=mm*f(1:9,i,j);
        %计算平衡动量
        m_eq(9,1)=rho(i,j);
        m_eq(1,1)=rho(i,j)*(-2+3*(u(i,j)^2+v(i,j)^2));
        m_eq(2,1)=rho(i,j)*(1-3*(u(i,j)^2+v(i,j)^2));
        m_eq(3,1)=rho(i,j)*u(i,j);
        m_eq(4,1)=-rho(i,j)*u(i,j);
        m_eq(5,1)=rho(i,j)*v(i,j);
        m_eq(6,1)=-rho(i,j)*v(i,j);
        m_eq(7,1)=rho(i,j)*(u(i,j)^2-v(i,j)^2);
        m_eq(8,1)=rho(i,j)*u(i,j)*v(i,j);
        %施加外力，需处理u计算，guo书；
        Fi=zeros(9,nx,ny);
        Fi(9,i,j)=0;
        Fi(1,i,j)=6*(1-0.5*s1(1))*(u(i,j)*forcex+v(i,j)*forcey);
        Fi(2,i,j)=-6*(1-0.5*s1(2))*(u(i,j)*forcex+v(i,j)*forcey);
        Fi(3,i,j)=(1-0.5*s1(3))*forcex;
        Fi(4,i,j)=-(1-0.5*s1(4))*forcex;
        Fi(5,i,j)=(1-0.5*s1(5))*forcey;
        Fi(6,i,j)=-(1-0.5*s1(6))*forcey;
        Fi(7,i,j)=2*(1-0.5*s1(7))*(u(i,j)*forcex-v(i,j)*forcey);
        Fi(8,i,j)=(1-0.5*s1(8))*(u(i,j)*forcey+v(i,j)*forcex);
        %在动量空间松弛
        m_temp=msl*(m-m_eq)-mminv*Fi(:,i,j);
        for k=1:9
            f(k,i,j)=f(k,i,j)-m_temp(k);
        end
    end
end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function result(nx,ny,u,v,count)
pi=0;
xn=[1:nx]'-1; yn=[1:ny]'-1;
for i=1:nx
    for j=1:ny
        pi=pi+1;
        peess(pi,:)=[xn(i),yn(j),u(i,j),v(i,j)];
    end
end
filename=['F:\LBM_code\date-1\' num2str(count) '-tecplot2d.dat'];
%   address是储存位置，这里的num2str是为了在循环输出dat数据文件中使用，如果只有一个文件可以忽略
fid=fopen(filename,'wt');
fprintf(fid,'variables= "x", "y", "U", "V"\r\n');
fprintf(fid,'zone t="Frame 0"i=%d,j=%d,f=point\r\n',ny,nx);
fprintf(fid,'SOLUTIONTIME=%d\r\n',count);
fprintf(fid,'%8.4f %8.4f %8.4f %8.4f \r\n',peess');
fclose(fid);
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function[rho,u,v]=ruv(f,nx,ny,cx,cy,forcex,forcey)
%calculate velocity compnents
for i=1:nx
    for j=1:ny
        rho(i,j)=sum(f(:,i,j));
        Usum=0.0; Vsum=0.0;
        for k=1:9
            Usum=Usum+f(k,i,j)*cx(k);
            Vsum=Vsum+f(k,i,j)*cy(k);
        end
        u(i,j)=Usum/rho(i,j)+forcex*0.5;
        v(i,j)=Vsum/rho(i,j)+forcey*0.5;
    end
end
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Streaming:
function [f]=stream(f,cx,cy)
for k=1:9
    f(k,:,:)=circshift( f(k,:,:), [0,cx(k),cy(k)] );
end
end
% End of streaming