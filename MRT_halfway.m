%LBM-MRT 2-D2Q9, note that c2=1/3, w9=4/9,
% 形式：f=zeros(nx,ny,9);
% w1-4=1/9, and w5-w8, 1/36
clear
nx=20;ny=50;
f=zeros(nx,ny,9);feq=zeros(nx,ny,9);fpre=zeros(nx,ny,9);
u=zeros(nx,ny);v=zeros(nx,ny);
rho=ones(nx,ny);  obst=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
c2=1./3.;
dx=1.0;dy=1.0;
x=(0:nx);
y=(0:ny);
L=ny;
% x0=nx/2+1; x1=ny/2+1;
% r0=nx/5;
% for i=1:nx
%     for j=1:ny
%         if (i-x0)^2+(j-x1)^2<=r0^2
%             obst(i,j)=1;
%         end
%     end
% end
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
             f(i,j,k)=rho(i,j)*w(k);%  计算初始态下的  分布函数=平衡分布函数                        
             Fi(k,i,j)=0;   %初始化离散的边界作用力分布函数,初始值赋为 0
         end
      end
 end
%Main Loop
tic
while error>tol
    % Collitions
    [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,s1,forcex,forcey);
    fpre=f;
    % Streaming:
    [f]=stream(f,cx,cy);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,fpre);
    
    % Calculate rho, u, v
    [rho,u,v]=ruv(f,forcex,forcey);
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
function [f]=boundary(nx,ny,f,fpre)
%bottom, and top boundary, bounce back
for i=1:nx
    f(i,1,2)=fpre(i,1,4);
    f(i,1,5)=fpre(i,1,7);
    f(i,1,6)=fpre(i,1,8);
    f(i,ny,4)=fpre(i,ny,2);
    f(i,ny,7)=fpre(i,ny,5);
    f(i,ny,8)=fpre(i,ny,6);
end
%Left boundary, Periodic
for j=1:ny
    f(1,j,1)=f(nx-1,j,1);
    f(1,j,5)=f(nx-1,j,5);
    f(1,j,8)=f(nx-1,j,8);
    f(nx,j,3)=f(2,j,3);
    f(nx,j,7)=f(2,j,7);
    f(nx,j,6)=f(2,j,6);
end

% for j=2:ny-1
%     f(nx,j,3)=f(nx-1,j,3);
%     f(nx,j,6)=f(nx-1,j,6);
%     f(nx,j,7)=f(nx-1,j,7);
%     f(1,j,1)=f(1,j,3)+2.*rho(1,j)*uo/3.;
%     f(1,j,5)=f(1,j,7)-0.5*(f(1,j,2)-f(1,j,4))+rho(1,j)*uo/6.;
%     f(1,j,8)=f(1,j,6)+0.5*(f(1,j,2)-f(1,j,4))+rho(1,j)*uo/6.;
%     u(1,j)=uo; v(1,j)=0.0;
% end

end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Collition
function [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,s1,forcex,forcey)
% Calculate the
f2=zeros(9,nx,ny);
for i=1:nx
    for j=1:ny
        for k=1:9
            f2(k,i,j)=f(i,j,k);
        end
        m=mm*f2(:,i,j);
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
            f2(k,i,j)=f2(k,i,j)-m_temp(k);
            f(i,j,k)=f2(k,i,j);
        end
    end
end
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function result(nx,ny,u,v,count)
pi=0;
xn=[1:nx]'-0.5; yn=[1:ny]'-0.5;
for i=1:nx
    for j=1:ny
        pi=pi+1;
        peess(pi,:)=[xn(i),yn(j),u(i,j),v(i,j)];
    end
end
filename=['F:\date\' num2str(count) '-tecplot2d.dat'];
%   address是储存位置，这里的num2str是为了在循环输出dat数据文件中使用，如果只有一个文件可以忽略
fid=fopen(filename,'wt');
fprintf(fid,'variables= "x", "y", "U", "V"\r\n');
fprintf(fid,'zone t="Frame 0"i=%d,j=%d,f=point\r\n',ny,nx);
fprintf(fid,'SOLUTIONTIME=%d\r\n',count);
fprintf(fid,'%8.4f %8.4f %8.4f %8.4f \r\n',peess');
fclose(fid);
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function[rho,u,v]=ruv(f,forcex,forcey)
rho=sum (f,3);
%calculate velocity compnents
u = ( sum(f(:,:,[1 5 8]),3) - sum(f(:,:,[3 6 7]),3) )./rho+forcex*0.5;
v = ( sum(f(:,:,[2 5 6]),3) - sum(f(:,:,[4 7 8]),3) )./rho+forcey*0.5;
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Streaming:
function [f]=stream(f,cx,cy)
for k=1:9
    f(:,:,k)=circshift(f(:,:,k), [cx(k),cy(k),0] );
end
end
% End of streaming