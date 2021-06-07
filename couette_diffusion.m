clear all; clc;
nx=31;ny=33;
Maxcount=400000;
f=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=zeros(nx,ny);
u=zeros(nx,ny);v=zeros(nx,ny);obst=zeros(nx,ny);
rho=ones(nx,ny);Cg=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];
c2=1./3.;
dx=1.0;dy=1.0;
x1=(0:1:nx); y1=(0:1:ny);
uo=0.10;
vo=0.10;
Re=10;
alpha=uo*(ny-1)/Re; 
omega=1./(3.*alpha+0.5);
Ma=1/(ny-1)/(3^0.5)*(1/omega-0.5)*Re;  
Pe=10;  
alphag=uo*(ny-1)/Pe; 
omegag=1.0./(3.*alphag+0.5);  
forcex=0; % 8.*uo*alpha./((ny-1).^2);
forcey=0;
%inital velocity
u(:,ny)=uo;
v(:,ny)=vo;
v(:,1)=vo;
count=0; cycle=40;
tol=1.0e-8; error=10.;erso=0.0;
tic;
while error>tol
    [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w,forcex,forcey);
    [f]=stream(f);
    [f]=boundary(nx,ny,f,uo,rho,u,v);
    [f]=obstc(nx,ny,f,opp,obst,u,v);
    [rho,u,v]=ruv(nx,ny,f);
    [g]=gcol(nx,ny,u,v,cx,cy,omegag,g,rhog,w); 
    [g]=stream(g); 
    [g]=gbound(nx,ny,w,g); 
    [g]=gobstc(nx,ny,g,opp,obst);
    rhog=sum(g,3);
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j)+rhog(i,j)*rhog(i,j);
        end
    end
    error=abs(ers-erso)/ers;
    erso=ers;
    if mod(count,cycle)==0
        result(nx,ny,u,v,rhog,count,obst);
    end
    count=count+1;
    if count>Maxcount
        break
    end
end
toc;
function [f]=boundary(nx,ny,f,uo,rho,u,v)
for j=2:ny-1
    f(nx,j,3)=f(1,j,3);
    f(nx,j,7)=f(1,j,7);
    f(nx,j,6)=f(1,j,6);
    f(1,j,1)=f(nx,j,1);
    f(1,j,5)=f(nx,j,5);
    f(1,j,8)=f(nx,j,8);
end
uo=0.1;
vo=0.1;
for i=1:nx
    %Top
    rhon=1./(1+vo)*(f(i,ny,9)+f(i,ny,1)+f(i,ny,3)+2.*(f(i,ny,2)+f(i,ny,6)+f(i,ny,5)));
    f(i,ny,4)=f(i,ny,2)-2./3*rhon*vo;
    f(i,ny,8)=f(i,ny,6)+1./2*(f(i,ny,3)-f(i,ny,1))+rhon*uo/2.0-1./6*rhon*vo;
    f(i,ny,7)=f(i,ny,5)+1./2*(f(i,ny,1)-f(i,ny,3))-rhon*uo/2.0-1./6*rhon*vo;
    %Below
    rhos=1./(1.-vo)*(f(i,1,9)+f(i,1,1)+f(i,1,3)+2.*(f(i,1,4)+f(i,1,8)+f(i,1,7)));
    f(i,1,2)=f(i,1,4)+2./3*rhos*vo;
    f(i,1,6)=f(i,1,8)+1./2*(f(i,1,1)-f(i,1,3))+rhos*vo/6.0;
    f(i,1,5)=f(i,1,7)+1./2*(f(i,1,3)-f(i,1,1))+rhon*uo/6.0;
end
end
function [f]=collition(nx,ny,u,v,cx,cy,omega,f,rho,w,forcex,forcey)
feq=zeros(nx,ny);
force=zeros(9);
for j=1:ny
    for i=1:nx
        t1=u(i,j)*u(i,j)+v(i,j)*v(i,j);
        for k=1:9
            t2=u(i,j)*cx(k)+v(i,j)*cy(k);
            force(k)=w(k).*3.*(cx(k).*forcex+cy(k).*forcey);
            feq(i,j,k)=rho(i,j)*w(k)*(1.0+3.0*t2+4.5*t2*t2-1.5*t1);
            f(i,j,k)=(1.-omega)*f(i,j,k)+omega*feq(i,j,k)+force(k);
        end
    end
end
end
function [g]=gbound(nx,ny,w,g)
g(1,:,1)=g(nx,:,1);
g(1,:,5)=g(nx,:,5);
g(1,:,8)=g(nx,:,8);
g(nx,:,3)=g(1,:,3);
g(nx,:,7)=g(1,:,7);
g(nx,:,6)=g(1,:,6);
co=1.0;
g(:,1,2)=(w(2)+w(4))*co-g(:,1,4);
g(:,1,5)=(w(5)+w(7))*co-g(:,1,7);
g(:,1,6)=(w(6)+w(8))*co-g(:,1,8);
c1=1.1;
g(:,ny,4)=(w(2)+w(4))*c1-g(:,ny,2);
g(:,ny,8)=(w(6)+w(8))*c1-g(:,ny,6);
g(:,ny,7)=(w(7)+w(5))*c1-g(:,ny,5);
end
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
filename=['F:\LBM_code\date-22\' num2str(count) '-tecplot2d.dat'];
fid=fopen(filename,'wt');
fprintf(fid,'variables= "x", "y", "U", "V", "C"\r\n');
fprintf(fid,'zone t="Frame 0"i=%d,j=%d,f=point\r\n',ny,nx);
fprintf(fid,'SOLUTIONTIME=%d\r\n',count);
fprintf(fid,'%8.4f %8.4f %8.4f %8.4f %8.4f \r\n',peess');
fclose(fid);
end
function[rho,u,v]=ruv(nx,ny,f)
rho=sum (f,3);
u = ( sum(f(:,:,[1 5 8]),3) - sum(f(:,:,[3 6 7]),3) )./rho;
v = ( sum(f(:,:,[2 5 6]),3) - sum(f(:,:,[4 7 8]),3) )./rho;
end
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
