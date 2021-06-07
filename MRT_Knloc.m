%LBM-MRT 2-D2Q9, note that c2=1/3, w9=4/9, ʧ��
% ��ʽ��f=zeros(nx,ny,9);
% w1-4=1/9, and w5-w8, 1/36
clear
nx=10;ny=10;
f=zeros(nx,ny,9);feq=zeros(nx,ny,9);fpre=zeros(nx,ny,9);
u=zeros(nx,ny);v=zeros(nx,ny);
Kn=zeros(nx,ny); sl=cell(nx,ny); Loc=zeros(nx,ny);obst=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
c2=1./3.;  arf=1.0;
L=ny;
Kn_out=0.0194; rho_out=3.00; rho_in=3.03;
rho=(rho_out+rho_in)/2*ones(nx,ny);
B1=1-0.1817*arf;   B2=0.8;
forcex=0.0;
forcey=0.0;
% uo=forcex*L^2/8/alpha2;
% forcey=0.0;
% Re=uo*L/alpha2;
% Ma=1/L/(3^0.5)*(ts2-0.5)*Re 
[Loc] = Llocal(obst,nx,ny,Loc);
M=[1 1 1 1 1 1 1 1 1;-4 -1 -1 -1 -1 2 2 2 2;4 -2 -2 -2 -2 1 1 1 1;...
    0 1 0 -1 0 1 -1 -1 1;0 -2 0 2 0 1 -1 -1 1;0 0 1 0 -1 1 1 -1 -1;...
    0 0 -2 0 2 1 1 -1 -1;0 1 -1 1 -1 0 0 0 0;0 0 0 0 0 1 -1 1 -1;];
mm=circshift(M,[-1 -1]);
mminv=inv(mm);
count=0; cycle=50;
tol=1.0e-6; error=10.;erso=0.0;
%��ʼ��
for i=1:nx
     for j=1:ny
          for k=1:9
             f(i,j,k)=rho(i,j)*w(k);%  �����ʼ̬�µ�  �ֲ�����=ƽ��ֲ�����                        
             Fi(k,i,j)=0;   %��ʼ����ɢ�ı߽��������ֲ�����,��ʼֵ��Ϊ 0
         end
      end
 end
%Main Loop
tic
while error>tol
    for i=1:nx
        for j=1:ny
            Kn(i,j)=Kn_out*rho_out*L/rho(i,j)/Loc(i,j);
            ts(i,j)=L*(6/pi)^0.5*Kn(i,j)./(1+2*Kn(i,j))+0.5;
            tq(i,j)=0.5+(3+4*pi*(ts(i,j)-0.5)^2*B2)/16./(ts(i,j)-0.5);
            omega(i,j)=1/ts(i,j);
            omegaq(i,j)=1/tq(i,j);
            r(i,j)=1/(1+B1*arf*(pi/6)^0.5);
            alpha(i,j)=c2*(1./omega(i,j)-0.5);
            sl{i,j}=[1.1 1.2 1 omegaq(i,j) 1 omegaq(i,j) omega(i,j) omega(i,j) 1];
            msl{i,j}=mminv*diag(sl{i,j});
        end
    end
    % Collitions
    [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,sl,forcex,forcey,Fi);
    fpre=f;
    % Streaming:
    [f]=stream(f,cx,cy);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,fpre,r,rho,rho_in,rho_out);
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
toc;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [f]=boundary(nx,ny,f,fpre,r,rho,rho_in,rho_out)
% %bottom, and top boundary, bounce back
for i=1:nx
    j=1;
    f(i,1,2)=f(i,1,4);
    f(i,1,5)=r(i,j)*fpre(i,1,7)+(1-r(i,j))*fpre(i,1,8);
    f(i,1,6)=r(i,j)*fpre(i,1,8)+(1-r(i,j))*fpre(i,1,7);
    j=ny;
    f(i,ny,4)=f(i,ny,2);
    f(i,ny,7)=r(i,j)*fpre(i,ny,5)+(1-r(i,j))*fpre(i,ny,6);
    f(i,ny,8)=r(i,j)*fpre(i,ny,6)+(1-r(i,j))*fpre(i,ny,5);
end
rho(1,:)=rho_in;
rho(nx,:)=rho_out;
for j=1:ny
    u(1,j)=1-(f(1,j,9)+f(1,j,2)+f(1,j,4)+2*(f(1,j,3)+f(1,j,6)+f(1,j,7)))/rho(1,j);
    f(1,j,1)=f(1,j,3)+2/3*rho(1,j)*u(1,j);
    f(1,j,5)=f(1,j,7)-0.5*(f(1,j,2)-f(1,j,4))+1/6*rho(1,j)*u(1,j);
    f(1,j,8)=f(1,j,6)+0.5*(f(1,j,2)-f(1,j,4))+1/6*rho(1,j)*u(1,j);

    u(nx,j)=(f(nx,j,9)+f(nx,j,2)+f(nx,j,4)+2*(f(nx,j,1)+f(nx,j,5)+f(nx,j,8)))/rho(nx,j)-1;
    f(nx,j,3)=f(nx,j,1)-2/3*rho(nx,j)*u(nx,j);
    f(nx,j,7)=f(nx,j,5)+0.5*(f(nx,j,2)-f(nx,j,4))-1/6*rho(nx,j)*u(nx,j);
    f(nx,j,6)=f(nx,j,8)-0.5*(f(nx,j,2)-f(nx,j,4))-1/6*rho(nx,j)*u(nx,j);
end
end
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Collition
function [f]=collition(nx,ny,u,v,f,rho,mm,mminv,msl,sl,forcex,forcey,Fi)
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
%����ƽ�⶯��
m_eq(9,:,:)=rho;
m_eq(1,:,:)=rho.*(-2+3.*(u.^2+v.^2));
m_eq(2,:,:)=rho.*(1-3.*(u.^2+v.^2));
m_eq(3,:,:)=rho.*u;
m_eq(4,:,:)=-rho.*u;
m_eq(5,:,:)=rho.*v;
m_eq(6,:,:)=-rho.*v;
m_eq(7,:,:)=rho.*(u.^2-v.^2);
m_eq(8,:,:)=rho.*u.*v;
for i=1:nx
    for j=1:ny
        Fi(9,i,j)=0;
        Fi(1,i,j)=6*(1-0.5*sl{i,j}(1))*(u(i,j)*forcex+v(i,j)*forcey);
        Fi(2,i,j)=-6*(1-0.5*sl{i,j}(2))*(u(i,j)*forcex+v(i,j)*forcey);
        Fi(3,i,j)=(1-0.5*sl{i,j}(3))*forcex;
        Fi(4,i,j)=-(1-0.5*sl{i,j}(4))*forcex;
        Fi(5,i,j)=(1-0.5*sl{i,j}(5))*forcey;
        Fi(6,i,j)=-(1-0.5*sl{i,j}(6))*forcey;
        Fi(7,i,j)=2*(1-0.5*sl{i,j}(7))*(u(i,j)*forcex-v(i,j)*forcey);
        Fi(8,i,j)=(1-0.5*sl{i,j}(8))*(u(i,j)*forcey+v(i,j)*forcex);
        m_temp{i,j}=msl{i,j}*(m(1:9,:)-m_eq(1:9,:))-mminv*Fi(1:9,:);
        mte=reshape(m_temp{i,j},9,nx,ny);
        for k=1:9
            f2(k,i,j)=f2(k,i,j)-mte(k,i,j);
        end
        
    end
end
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
%   address�Ǵ���λ�ã������num2str��Ϊ����ѭ�����dat�����ļ���ʹ�ã����ֻ��һ���ļ����Ժ���
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