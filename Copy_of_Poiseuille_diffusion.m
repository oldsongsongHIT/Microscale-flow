%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear all; clc;
nx=5;ny=32;
Maxcount=400000;
f=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=zeros(nx,ny);
u=zeros(nx,ny);v=zeros(nx,ny);obst=zeros(nx,ny);
rho=ones(nx,ny);Cg=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];
c2=1./3.;
uo=0.10; 
Kn=0.1;
omega=1./(Kn*(ny-1)*(6/pi)^0.5+0.5);
alpha=c2*(1/omega-0.5);
r=1;  %���Ʊ߽�rfϵ��
pe=10;  %佩克莱特�?
alphag=uo*(ny-1)/pe; %扩散系数
omegag=1.0./(3.*alphag+0.5);  %扩散系数对应松弛时间
%force
forcex=8.*uo*alpha./((ny-1).^2);
forcey=0;
count=0; cycle=40;%每隔cycle输出�?次结�?
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
% addpath(genpath('E:\LBM\Code\Reaction_flow'));
tic;
while error>tol
    % Collitions
    [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w,forcex,forcey);
    % Streaming:
    [f]=stream(f);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,uo,rho,u,v,r);
    %Obsticale-流动
     [f]=obstc(nx,ny,f,opp,obst,u,v);
    % Calculate rho, u, v
    [rho,u,v]=ruv(nx,ny,f,forcex,forcey);
%     %溶质组分
%     [g]=gcol(nx,ny,u,v,cx,cy,omegag,g,rhog,w); %collision
%     [g]=stream(g); %streaming
%     [g]=gbound(nx,ny,w,g); % boundary conditions
%     %组分扩散-Obstacle
%     [g]=gobstc(nx,ny,g,opp,obst);
%     rhog=sum(g,3);%%计算浓度
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j)+rhog(i,j)*rhog(i,j);
        end
    end
    error=abs(ers-erso)/ers;
    erso=ers;
    %输出文件
    if mod(count,cycle)==0
        result(nx,ny,u,v,rhog,count,obst);
    end
    count=count+1;
    if count>Maxcount
        break
    end
end
Uz=u./uo;
toc;
%各对应函�?
% %Boudary conditions
function [f]=boundary(nx,ny,f,uo,rho,u,v,r)
%right hand boundary, Periodic
for j=1:ny
    f(nx,j,3)=f(1,j,3);
    f(nx,j,7)=f(1,j,7);
    f(nx,j,6)=f(1,j,6);
end
%bottom, and top boundary, bounce back
% for i=1:nx
%     f(i,1,2)=f(i,1,4);
%     f(i,1,5)=f(i,1,7);
%     f(i,1,6)=f(i,1,8);
%     f(i,ny,4)=f(i,ny,2);
%     f(i,ny,7)=f(i,ny,5);
%     f(i,ny,8)=f(i,ny,6);
% end

f(:,1,2)=f(:,1,4);
f(:,1,5)=r*f(:,1,7)+(1-r)*f(:,1,8);
f(:,1,6)=r*f(:,1,8)+(1-r)*f(:,1,7);
f(:,ny,4)=f(:,ny,2);
f(:,ny,7)=r*f(:,ny,5)+(1-r)*f(:,ny,8);
f(:,ny,8)=r*f(:,ny,6)+(1-r)*f(:,ny,7);


%Left boundary, Periodic
for j=2:ny-1
    f(1,j,1)=f(nx,j,1);
    f(1,j,5)=f(nx,j,5);
    f(1,j,8)=f(nx,j,8);
end
% End of boundary conditions.
end
% Collition
function [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w,forcex,forcey)
feq=zeros(nx,ny);
force=zeros(9);
for j=1:ny
    for i=1:nx
        t1=u(i,j)*u(i,j)+v(i,j)*v(i,j);
        for k=1:9
            t2=u(i,j)*cx(k)+v(i,j)*cy(k);
%             force(k)=3.*w(k)*(cx(k).*forcex+cy(k).*forcey);
%              force(k)=(1.-0.5*omega)*w(k)*rho(i,j).*...
%              (3.*((cx(k)-u(i,j)).*forcex+(cy(k)-v(i,j)).*forcey)+...
%              9*(t2*(cx(k)*forcex+cy(k)*forcey)));
            feq(i,j,k)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
            force(k)=feq(i,j,k)*3*((cx(k)-u(i,j))*forcex+(cy(k)-v(i,j))*forcey);
            f(i,j,k)=(1.-omega)*f(i,j,k)+omega*feq(i,j,k)+force(k);
        end
    end
end
end
%组分对流方程的边界条�?
function [g]=gbound(nx,ny,w,g)
%Boundary condition:
%left boundary, iPeriodic
g(1,:,1)=g(nx,:,1);
g(1,:,5)=g(nx,:,5);
g(1,:,8)=g(nx,:,8);
%right hand boundary,  Periodic
g(nx,:,3)=g(1,:,3);
g(nx,:,7)=g(1,:,7);
g(nx,:,6)=g(1,:,6);
%bottom boundary, co
co=1.0;
g(:,1,2)=(w(2)+w(4))*co-g(:,1,4);
g(:,1,5)=(w(5)+w(7))*co-g(:,1,7);
g(:,1,6)=(w(6)+w(8))*co-g(:,1,8);
%Top boundary, c1
c1=1.1;
g(:,ny,4)=(w(2)+w(4))*c1-g(:,ny,2);
g(:,ny,8)=(w(6)+w(8))*c1-g(:,ny,6);
g(:,ny,7)=(w(7)+w(5))*c1-g(:,ny,5);
end
%Obsticale replace at the entrance
function [f]=obstc(nx,ny,f,opp,obst,u,v)
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            for k=1:9
                f(i,j,k) = f(i,j,opp(k));
            end
            u(i,j)=0.0;
            v(i,j)=0.0;
        end
    end
end
end
% 组分对流方程的Collition
function [g]=gcol(nx,ny,u,v,cx,cy,omegag,g,rhog,w)
geq=zeros(nx,ny);
for j=1:ny
    for i=1:nx
        for k=1:9
            t=u(i,j)*cx(k)+v(i,j)*cy(k);
            geq(i,j,k)=rhog(i,j)*w(k)*(1.0+3.0*t);
            g(i,j,k)=(1.-omegag)*g(i,j,k)+omegag*geq(i,j,k);
        end
    end
end
end
% Obsticale-组分扩散方程障碍
function [g]=gobstc(nx,ny,g,opp,obst)
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            for k=1:9
                g(i,j,k) = g(i,j,opp(k));
            end
        end
    end
end
end
function result(nx,ny,u,v,rhog,count,obst)
Cg=rhog;
pi=0;
xn=[1:nx]'; yn=[1:ny]';
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            Cg(i,j)=0.0;
        end
        pi=pi+1;
        peess(pi,:)=[xn(i),yn(j),u(i,j),v(i,j),Cg(i,j)];
    end
end
filename=['F:\LBM_code\date-2\' num2str(count) '-tecplot2d.dat'];
%   address是储存位置，这里的num2str是为了在循环输出dat数据文件中使用，如果只有�?个文件可以忽�?
fid=fopen(filename,'wt');
fprintf(fid,'variables= "x", "y", "U", "V", "C"\r\n');
fprintf(fid,'zone t="Frame 0"i=%d,j=%d,f=point\r\n',ny,nx);
fprintf(fid,'SOLUTIONTIME=%d\r\n',count);
fprintf(fid,'%8.4f %8.4f %8.4f %8.4f %8.4f \r\n',peess');
fclose(fid);
end
function[rho,u,v]=ruv(nx,ny,f,forcex,forcey)
rho=sum (f,3);
%calculate velocity compnents
% u = ( sum(f(:,:,[1 5 8]),3) - sum(f(:,:,[3 6 7]),3) +forcex*2./2)./rho;
% v = ( sum(f(:,:,[2 5 6]),3) - sum(f(:,:,[4 7 8]),3) +forcey*2./2)./rho;


u = ( sum(f(:,:,[1 5 8]),3) - sum(f(:,:,[3 6 7]),3) )./rho;
v = ( sum(f(:,:,[2 5 6]),3) - sum(f(:,:,[4 7 8]),3) )./rho;
end
% Streaming:
function [f]=stream(f)
f(:,:,1)=circshift( squeeze(f(:,:,1)), [+1,+0] );
f(:,:,2)=circshift( squeeze(f(:,:,2)), [+0,+1] );
f(:,:,3)=circshift( squeeze(f(:,:,3)), [-1,+0] );
f(:,:,4)=circshift( squeeze(f(:,:,4)), [+0,-1] );
f(:,:,5)=circshift( squeeze(f(:,:,5)), [+1,+1] );
f(:,:,6)=circshift( squeeze(f(:,:,6)), [-1,+1] );
f(:,:,7)=circshift( squeeze(f(:,:,7)), [-1,-1] );
f(:,:,8)=circshift( squeeze(f(:,:,8)), [+1,-1] );
end
% End of streaming
