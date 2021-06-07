%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear;
nx=33;ny=65;
Maxcount=400000;
f=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=zeros(nx,ny);feq=zeros(nx,ny,9);
u=zeros(nx,ny);v=zeros(nx,ny);obst=zeros(nx,ny);
rho=3*ones(nx,ny);Kn=zeros(nx,ny); omega=zeros(nx,ny);r=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];
c2=1./3.;
dx=1.0;dy=1.0;
x1=(0:1:nx); y1=(0:1:ny);
L=ny-1;
uo=0.00; 
pe=10;  %佩克莱特�?
alphag=uo*L/pe; %扩散系数
omegag=1.0./(3.*alphag+0.5);  %扩散系数对应松弛时间
%force
forcex=0.0;
forcey=0;
count=0; cycle=40;
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
tic;
while error>tol
    Kn_out=0.194;  rho_out=3.0;
    Kn=Kn_out*rho_out./rho;
    tf=L*(6/pi)^0.5.*Kn+0.5;
    teff=0.5+(tf-0.5).*(1./(1+2.2*Kn).*rho/rho_out);
    %(2/pi.*atan(2^0.5.*Kn.^(-0.75)));
    omega=1./teff;
    alpha=c2*(1./omega-0.5);
    arf=1;
    A1=(2-arf)/arf*(1-0.1817*arf);
    A2=arf^2*(1/pi+0.5*A1^2);
    r=1./(1+(pi/6)^0.5.*(1/4./Kn./L^2.+A1+(2*A2*Kn-8./pi.*Kn)));  %���Ʊ߽�rfϵ��
    % Collitions
    [f]=collition(nx,ny,u,v,cx,cy,omega,f,feq,rho,w,forcex,forcey);
    % Streaming:
    [f]=stream(f);
    % End of streaming
    %Boundary condition:
    [f,u,v,rho]=boundary(nx,ny,f,uo,rho,u,v,r);
    % Calculate rho, u, v
    [rho,u,v]=ruv(nx,ny,f,forcex,forcey);
    %溶质组分
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j);
        end
    end
    error=abs(ers-erso)/ers;
    erso=ers;
    %输出文件
    if mod(count,cycle)==0
        result(nx,ny,u,v,rhog,count,obst,uo);
    end
    count=count+1;
    if count>Maxcount
        break
    end
end
toc;
%
function [f]=collition(nx,ny,u,v,cx,cy,omega,f,feq,rho,w,forcex,forcey)
force=zeros(9);
for j=1:ny
    for i=1:nx
        t1=u(i,j)*u(i,j)+v(i,j)*v(i,j);
        for k=1:9
            t2=u(i,j)*cx(k)+v(i,j)*cy(k);
             force(k)=(1.-0.5*omega(i,j))*w(k)*rho(i,j).*(3.*((cx(k)-u(i,j)).*forcex+(cy(k)-v(i,j)).*forcey)+9*(t2*(cx(k)*forcex+cy(k)*forcey)));
             feq(i,j,k)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
%             force(k)=feq(i,j,k)*3*((cx(k)-u(i,j))*forcex+(cy(k)-v(i,j))*forcey);
            f(i,j,k)=(1.-omega(i,j))*f(i,j,k)+omega(i,j)*feq(i,j,k)+force(k);
        end
    end
end
end
% %Boudary conditions
function [f,u,v,rho]=boundary(nx,ny,f,uo,rho,u,v,r)
%right hand boundary, Periodic
Pin=1.3;
rho(1,:)=3*Pin;
Pout=1.00;
rho(nx,:)=3*Pout;
for j=2:ny-1
    u(1,j)=1-(f(1,j,9)+f(1,j,2)+f(1,j,4)+2*(f(1,j,3)+f(1,j,6)+f(1,j,7)))/rho(1,j);
    f(1,j,1)=f(1,j,3)+2/3*rho(1,j)*u(1,j);
    f(1,j,5)=f(1,j,7)-0.5*(f(1,j,2)-f(1,j,4))+1/6*rho(1,j)*u(1,j);
    f(1,j,8)=f(1,j,6)+0.5*(f(1,j,2)-f(1,j,4))+1/6*rho(1,j)*u(1,j);

    u(nx,j)=(f(nx,j,9)+f(nx,j,2)+f(nx,j,4)+2*(f(nx,j,1)+f(nx,j,5)+f(nx,j,8)))/rho(nx,j)-1;
    f(nx,j,3)=f(nx,j,1)-2/3*rho(nx,j)*u(nx,j);
    f(nx,j,7)=f(nx,j,5)+0.5*(f(nx,j,2)-f(nx,j,4))-1/6*rho(nx,j)*u(nx,j);
    f(nx,j,6)=f(nx,j,8)-0.5*(f(nx,j,2)-f(nx,j,4))-1/6*rho(nx,j)*u(nx,j);
end
%bottom, and top boundary, bounce back
for i=2:nx-1
    j=1;
    f(i,1,2)=f(i,1,4);
    f(i,1,5)=r(i,j)*f(i,1,7)+(1-r(i,j))*f(i,1,8);
    f(i,1,6)=r(i,j)*f(i,1,8)+(1-r(i,j))*f(i,1,7);
    j=ny;
    f(i,ny,4)=f(i,ny,2);
    f(i,ny,7)=r(i,j)*f(i,ny,5)+(1-r(i,j))*f(i,ny,6);
    f(i,ny,8)=r(i,j)*f(i,ny,6)+(1-r(i,j))*f(i,ny,5);
end
% corner 1
f(1,1,1)=f(1,1,3);
f(1,1,2)=f(1,1,4);
f(1,1,5)=r(1,1)*f(1,1,7)+(1-r(1,1))*f(1,1,8);
% corner 2
f(1,ny,1)=f(1,ny,3);
f(1,ny,4)=f(1,ny,2);
f(1,ny,8)=r(1,ny)*f(1,ny,6)+(1-r(1,ny))*f(1,ny,5);
% corner 3
f(nx,1,3)=f(nx,1,1);
f(nx,1,2)=f(nx,1,4);
f(nx,1,6)=r(nx,1)*f(nx,1,8)+(1-r(nx,1))*f(nx,1,7);
% corner 4
f(nx,ny,3)=f(nx,ny,1);
f(nx,ny,4)=f(nx,ny,2);
f(nx,ny,7)=r(nx,ny)*f(nx,ny,5)+(1-r(nx,ny))*f(nx,ny,6);
% End of boundary conditions.
end

function result(nx,ny,u,v,rhog,count,obst,uo)
Uz=((u.^2+v.^2).^0.5)./uo;
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
