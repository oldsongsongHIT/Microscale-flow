%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear all;
nx=5;ny=501;
Maxcount=400000;
f=zeros(nx,ny,9);feq=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=zeros(nx,ny);
u=zeros(nx,ny);v=zeros(nx,ny);obst=zeros(nx,ny);
rho=ones(nx,ny);Cg=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];
c2=1./3.;
dx=1.0;dy=1.0;
x1=(0:1:nx); y1=(0:1:ny);
obst(:,[1,ny])=1.0;
L=ny-1;
arf=1;
Kn=0.1;
tf=(L*(6/pi)^0.5.*Kn+0.5);
omega=1./tf;
r=1;
% %guo £¿
% tf=(L*(6/pi)^0.5.*Kn+0.5).*(2/pi.*atan(2^0.5.*Kn.^(-0.75)));
% omega=1./tf;
% A1=(2-arf)/arf*(1-0.1817*arf);
% A2=arf^2*(1/pi+0.5*A1^2);
% r=1/(1+(pi/6)^0.5*(1/4/Kn/L^2+A1+(2*A2*Kn-8/pi*Kn)));  %»¬ÒÆ±ß½çrfÏµÊý
% % Liqing (2011)
% B1=1-0.1817*arf;   B2=0.8;
% tf=L*(6/pi)^0.5.*Kn/(1+2*Kn)+0.5;
% tq=0.5+(3+4*pi*(tf-0.5)^2*B2)/16/(tf-0.5);
% omega=1./tq;
% r=1/(1+B1*arf*(pi/6)^0.5);
%Zhang 2006
tf=Kn*L/(8/3/pi)^0.5+0.5;
%½áÊø
alpha=c2*(1/omega-0.5);
forcex=1e-4;
uo=forcex*L^2/8/alpha;
forcey=0.0;
Re=uo*L/alpha;
Ma=1/L/(3^0.5)*(1/omega-0.5)*Re
count=0; cycle=40;
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
tic;
while error>tol
    % Collitions
    [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w,forcex,forcey);
    % Streaming:
    [f]=stream(f);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,uo,rho,u,v,r);
%     [f]=obstc(nx,ny,f,uo,rho,opp,obst,r);
    % Calculate rho, u, v
    [rho,u,v]=ruv(nx,ny,f,forcex,forcey);
    %æº¶è´¨ç»„åˆ†
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j);
        end
    end
    error=abs(ers-erso)/ers;
    erso=ers;
    %è¾“å‡ºæ–‡ä»¶
%     if mod(count,cycle)==0
%         result(nx,ny,u,v,rhog,count,obst,uo);
%     end
    count=count+1;
    if count>Maxcount
        break
    end
end
Uz=u./uo;
toc;
%
function [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w,forcex,forcey)
feq=zeros(nx,ny);
force=zeros(9);
for j=1:ny
    for i=1:nx
        t1=u(i,j)*u(i,j)+v(i,j)*v(i,j);
        for k=1:9
            t2=u(i,j)*cx(k)+v(i,j)*cy(k);
             force(k)=(1.-0.5*omega)*w(k)*rho(i,j).*...
             (3.*((cx(k)-u(i,j)).*forcex+(cy(k)-v(i,j)).*forcey)+...
             9*(t2*(cx(k)*forcex+cy(k)*forcey)));
             feq(i,j,k)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
%             force(k)=feq(i,j,k)*3*((cx(k)-u(i,j))*forcex+(cy(k)-v(i,j))*forcey);
            f(i,j,k)=(1.-omega)*f(i,j,k)+omega*feq(i,j,k)+force(k);
        end
    end
end
end
% %Boudary conditions
function [f]=boundary(nx,ny,f,uo,rho,u,v,r)
%right hand boundary, Periodic
for j=1:ny
    f(nx,j,3)=f(2,j,3);
    f(nx,j,7)=f(2,j,7);
    f(nx,j,6)=f(2,j,6);
    f(1,j,1)=f(nx-1,j,1);
    f(1,j,5)=f(nx-1,j,5);
    f(1,j,8)=f(nx-1,j,8);
end
%bottom, and top boundary, bounce back
for i=1:nx
    f(i,1,2)=f(i,1,4);
    f(i,1,5)=r*f(i,1,7)+(1-r)*f(i,1,8);
    f(i,1,6)=r*f(i,1,8)+(1-r)*f(i,1,7);
    f(i,ny,4)=f(i,ny,2);
    f(i,ny,7)=r*f(i,ny,5)+(1-r)*f(i,ny,6);
    f(i,ny,8)=r*f(i,ny,6)+(1-r)*f(i,ny,5);
end
%Left boundary, Periodic

% End of boundary conditions.
end
% Collition

function [f]=obstc(nx,ny,f,uo,rho,opp,obst,r)
%length of obsticale= nx/5, and has sides of 10 units
obb=[8,7,6,5];
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            for k=1:4
                f(i,j,k) = f(i,j,opp(k));
            end
            for k=5:8
                f(i,j,k)=r*f(i,j,opp(k))+(1-r)*f(i,j,obb(k-4));
            end
 
        end
    end
end
end

function result(nx,ny,u,v,rhog,count,obst,uo)
pi=0;
xn=[1:nx]'; yn=[1:ny]';
for i=1:nx
    for j=1:ny
        pi=pi+1;
        peess(pi,:)=[xn(i),yn(j),Uz(i,j)];
    end
end
filename=['F:\LBM_code\date-1\' num2str(count) '-tecplot2d.dat'];
fid=fopen(filename,'wt');
fprintf(fid,'variables= "x", "y", "U"\n');
fprintf(fid,'zone t="Frame 0"i=%d,j=%d,f=point\n',ny,nx);
fprintf(fid,'SOLUTIONTIME=%d\n',count);
fprintf(fid,'%8.4f %8.4f %8.4f\n',peess');
fclose(fid);
end
function[rho,u,v]=ruv(nx,ny,f,forcex,forcey)
rho=sum (f,3);
%calculate velocity compnents
u = ( sum(f(:,:,[1 5 8]),3) - sum(f(:,:,[3 6 7]),3) +rho.*forcex*1./2)./rho;
v = ( sum(f(:,:,[2 5 6]),3) - sum(f(:,:,[4 7 8]),3) +rho.*forcey*1./2)./rho;
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
