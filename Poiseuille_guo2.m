%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear;
nx=5;ny=32;
Maxcount=400000;
f1=zeros(nx,ny,9);fpost=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=zeros(nx,ny);feq=zeros(nx,ny,9);
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
L=ny;
uo=0.10;  %�?般取0.1数量�?,�?般小�?0.3
Re=10;
alpha=uo*L/Re; %�?般取0.01数量�?,�?<0.1�?
omega=1./(3.*alpha+0.5);
Ma=1/L/(3^0.5)*(1/omega-0.5)*Re;  %�?般小�?0.2，越小越�?%
pe=10;  %佩克莱特�?
alphag=uo*L/pe; %扩散系数
omegag=1.0./(3.*alphag+0.5);  %扩散系数对应松弛时间
%force
forcex=8*uo*alpha/(L^2);
forcey=0;
count=0; cycle=40;
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
tic;
while error>tol
    % Collitions
    [fpost]=collition(nx,ny,u,v,cx,cy,omega,fpost,feq,rho,w,forcex,forcey);
    f1=fpost;
    % Streaming:
    [fpost]=stream(fpost);
    % End of streaming
    %Boundary condition:
    [fpost]=boundary(nx,ny,f1,fpost,uo,rho,u,v);
    
    [fpost]=obstc(nx,ny,fpost,f1,opp,obst,u,v);
    % Calculate rho, u, v
    [rho,u,v]=ruv(nx,ny,fpost,forcex,forcey);
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
Uz=u./uo;
toc;
%
function [fpost]=collition(nx,ny,u,v,cx,cy,omega,fpost,feq,rho,w,forcex,forcey)
force=zeros(9);
for j=1:ny
    for i=1:nx
        t1=u(i,j)*u(i,j)+v(i,j)*v(i,j);
        for k=1:9
            t2=u(i,j)*cx(k)+v(i,j)*cy(k);
             force(k)=(1.-0.5*omega)*w(k)*rho(i,j).*(3.*((cx(k)-u(i,j)).*forcex+(cy(k)-v(i,j)).*forcey)+9*(t2*(cx(k)*forcex+cy(k)*forcey)));
             feq(i,j,k)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
%             force(k)=feq(i,j,k)*3*((cx(k)-u(i,j))*forcex+(cy(k)-v(i,j))*forcey);
            fpost(i,j,k)=(1.-omega)*fpost(i,j,k)+omega*feq(i,j,k)+force(k);
        end
    end
end
end
% %Boudary conditions
function [fpost]=boundary(nx,ny,f1,fpost,uo,rho,u,v)
%right hand boundary, Periodic
for j=1:ny
    fpost(nx,j,3)=fpost(2,j,3);
    fpost(nx,j,7)=fpost(2,j,7);
    fpost(nx,j,6)=fpost(2,j,6);
end
%bottom, and top boundary, bounce back
% for i=1:nx
%     fpost(i,1,2)=f1(i,1,4);
%     fpost(i,1,5)=f1(i,1,7);
%     fpost(i,1,6)=f1(i,1,8);
%     fpost(i,ny,4)=f1(i,ny,2);
%     fpost(i,ny,7)=f1(i,ny,5);
%     fpost(i,ny,8)=f1(i,ny,6);
% end
%Left boundary, Periodic
for j=1:ny
    fpost(1,j,1)=fpost(nx-1,j,1);
    fpost(1,j,5)=fpost(nx-1,j,5);
    fpost(1,j,8)=fpost(nx-1,j,8);
end
% End of boundary conditions.
end

function [fpost]=obstc(nx,ny,fpost,f1,opp,obst,u,v)
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            for k=1:9
                fpost(i,j,k) = f1(i,j,opp(k));
            end
        end
    end
end
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
function[rho,u,v]=ruv(nx,ny,fpost,forcex,forcey)
rho=sum (fpost,3);
%calculate velocity compnents
u = ( sum(fpost(:,:,[1 5 8]),3) - sum(fpost(:,:,[3 6 7]),3) +rho.*forcex*1./2)./rho;
v = ( sum(fpost(:,:,[2 5 6]),3) - sum(fpost(:,:,[4 7 8]),3) +rho.*forcey*1./2)./rho;
end
% Streaming:
function [fpost]=stream(fpost)
fpost(:,:,1)=circshift( squeeze(fpost(:,:,1)), [+1,+0] );
fpost(:,:,2)=circshift( squeeze(fpost(:,:,2)), [+0,+1] );
fpost(:,:,3)=circshift( squeeze(fpost(:,:,3)), [-1,+0] );
fpost(:,:,4)=circshift( squeeze(fpost(:,:,4)), [+0,-1] );
fpost(:,:,5)=circshift( squeeze(fpost(:,:,5)), [+1,+1] );
fpost(:,:,6)=circshift( squeeze(fpost(:,:,6)), [-1,+1] );
fpost(:,:,7)=circshift( squeeze(fpost(:,:,7)), [-1,-1] );
fpost(:,:,8)=circshift( squeeze(fpost(:,:,8)), [+1,-1] );
end
% End of streaming
