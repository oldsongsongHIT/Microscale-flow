%LBM-MRT 2-D2Q9, note that c2=1/3, w9=4/9,
% 形式：f=zeros(nx,ny,9);
% w1-4=1/9, and w5-w8, 1/36
clear
nx=5;ny=50;
f=zeros(nx,ny,9);feq=zeros(nx,ny,9);fpre=zeros(nx,ny,9);
u=zeros(nx,ny);v=zeros(nx,ny);
rho=ones(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
c2=1./3.;  arf=1.0;
L=ny;
Kn_ref=0.1; 
lan_ref=Kn*L;
%积分
syms t;
E=t^(-1)*exp(-t);
y=L; %函数获得
a=y/lan_ref;
Ei=int(Ft,t,a,inf);
K=1+(a-1)*exp(-a)-a^2*Ei;
ts=L*(6/pi)^0.5.*Kn*K+0.5;
omega=1./ts;
r=zeros(nx,2); tq=zeros(nx,2);
A1=(2-arf)/arf*(1-0.1817*arf);
A2=1/pi+0.5*A1^2;
r(:,1)=1/(1+(pi/6)^0.5*A1+(ts(:,2)-ts(:,1))/8/(ts(:,1)-0.5));
tq(:,1)=0.5+(3+24*pi/6*(ts(:,1)-0.5)^2*A2)/16/(ts(:,1)-0.5)+((ts(:,2)-ts(:,1))*(12+30*(ts(:,1)-0.5)*(pi/6)^0.5*A1)/16/(ts(:,1)-0.5)^2);
r(:,ny)=1/(1+(pi/6)^0.5*A1+(ts(:,2)-ts(:,ny))/8/(ts(:,ny)-0.5));
tq(:,ny)=0.5+(3+24*pi/6*(ts(:,ny)-0.5)^2*A2)/16/(ts(:,ny)-0.5)+((ts(:,2)-ts(:,ny))*(12+30*(ts(:,ny)-0.5)*(pi/6)^0.5*A1)/16/(ts(:,ny)-0.5)^2);
omega=1./ts;
omegaq=1./tq;
%结束
alpha=c2*(1/omega-0.5);
forcex=1e-9;
uo=forcex*L^2/8/alpha;
forcey=0.0;
count=0; cycle=50;
Re=uo*L/alpha;
Ma=1/L/(3^0.5)*(1/omega-0.5)*Re 
tol=1.0e-20; error=10.;erso=0.0;
s1=[1.1 1.2 1 omegaq 1 omegaq omega omega 1];
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
    [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,s1,forcex,forcey,Fi);
    fpre=f;
    % Streaming:
    [f]=stream(f,cx,cy);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,fpre,r);
    % Calculate rho, u, v
    [rho,u,v]=ruv(f,forcex,forcey);
    count=count+1;
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j);
        end
    end
    error=abs(ers-erso)/ers;
    erso=ers;
    if mod(count,50)==0
        result(nx,ny,u,v,count);
    end
end
result(nx,ny,u,v,count);
Uz=u./uo;
toc;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [f]=boundary(nx,ny,f,fpre,r)
%bottom, and top boundary, bounce back
for i=1:nx
    f(i,1,2)=fpre(i,1,4);
    f(i,1,5)=r*fpre(i,1,7)+(1-r)*fpre(i,1,8);
    f(i,1,6)=r*fpre(i,1,8)+(1-r)*fpre(i,1,7);
    f(i,ny,4)=fpre(i,ny,2);
    f(i,ny,7)=r*fpre(i,ny,5)+(1-r)*fpre(i,ny,6);
    f(i,ny,8)=r*fpre(i,ny,6)+(1-r)*fpre(i,ny,5);
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
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Collition
function [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,s1,forcex,forcey,Fi)
% Calculate the
f2=zeros(9,nx,ny);
for i=1:nx
    for j=1:ny
        for k=1:9
            f2(k,i,j)=f(i,j,k);
        end
    end
end
m=mm*f2(1:9,:);
%计算平衡动量
m_eq(9,:,:)=rho;
m_eq(1,:,:)=rho.*(-2+3.*(u.^2+v.^2));
m_eq(2,:,:)=rho.*(1-3.*(u.^2+v.^2));
m_eq(3,:,:)=rho.*u;
m_eq(4,:,:)=-rho.*u;
m_eq(5,:,:)=rho.*v;
m_eq(6,:,:)=-rho.*v;
m_eq(7,:,:)=rho.*(u.^2-v.^2);
m_eq(8,:,:)=rho.*u.*v;
Fi(9,:,:)=0;
Fi(1,:,:)=6*(1-0.5*s1(1))*(u.*forcex+v.*forcey);
Fi(2,:,:)=-6*(1-0.5*s1(2))*(u.*forcex+v.*forcey);
Fi(3,i,j)=(1-0.5*s1(3))*forcex;
Fi(4,i,j)=-(1-0.5*s1(4))*forcex;
Fi(5,i,j)=(1-0.5*s1(5))*forcey;
Fi(6,:,:)=-(1-0.5*s1(6))*forcey;
Fi(7,:,:)=2*(1-0.5*s1(7))*(u.*forcex-v.*forcey);
Fi(8,:,:)=(1-0.5*s1(8))*(u.*forcey+v.*forcex);
%在动量空间松弛
m_temp=msl*(m(1:9,:)-m_eq(1:9,:))-mminv*Fi(1:9,:);
f2(:)=f2(:)-m_temp(:);
for i=1:nx
    for j=1:ny
        for k=1:9
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
filename=['F:\LBM_code\date-2\' num2str(count) '-tecplot2d.dat'];
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