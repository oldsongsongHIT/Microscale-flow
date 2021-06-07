%LBM-MRT 2-D2Q9, note that c2=1/3, w9=4/9,
% 形式：f=zeros(nx,ny,9);
% w1-4=1/9, and w5-w8, 1/36
clear
nx=5;ny=21;
f=zeros(nx,ny,9);feq=zeros(nx,ny,9);f=zeros(nx,ny,9);
u=zeros(nx,ny);v=zeros(nx,ny);
rho=ones(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
c2=1./3.;  arf=1;
L=ny-1;
Kn=0.5; 
% landa=Kn*L;
% ts=(L*(6/pi)^0.5.*Kn+0.5).*(2/pi.*atan(2^0.5.*Kn.^(-0.75)));
% omega=1./ts;
% A1=(2-arf)/arf*(1-0.1817*arf);
% A2=arf^2*(1/pi+0.5*A1^2);
% r=1/(1+(pi/6)^0.5*(1/4/Kn/L^2+A1+(2*A2*Kn-8/pi*Kn)));  %滑移边界rf系数
% % Liqing (2011)
B1=1-0.1817*arf;   B2=0.8;
ts=L*(6/pi)^0.5.*Kn/(1+2*Kn)+0.5;
tq=0.5+(3+4*pi*(ts-0.5)^2*B2)/16/(ts-0.5);
omega=1./ts;
omegaq=1./tq;
r=1/(1+B1*arf*(pi/6)^0.5);
%结束
alpha=c2*(1/omega-0.5);
forcex=1e-11;
uo=0.1;%forcex*L^2/8/alpha;
forcey=0.0;
count=0; cycle=40;
Re=uo*L/alpha;
Ma=1/L/(3^0.5)*(1/omega-0.5)*Re 
tol=1.0e-8; error=10.;erso=0.0;
s1=[omega omega-0.1 0 omegaq 0 omegaq omega omega 0];
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
    % Streaming:
    [f]=stream(f,cx,cy);
    % End of streaming
    %Boundary condition:
    [f]=boundary(nx,ny,f,r,rho,w,cx,cy,c2);
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
Uz=u./uo;
toc;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [f]=boundary(nx,ny,f,r,rho,w,cx,cy,c2)
%bottom, and top boundary, bounce back
for j=1:ny
    f(1,j,1)=f(nx-1,j,1);
    f(1,j,5)=f(nx-1,j,5);
    f(1,j,8)=f(nx-1,j,8);
    f(nx,j,3)=f(2,j,3);
    f(nx,j,7)=f(2,j,7);
    f(nx,j,6)=f(2,j,6);
end
uw=0.1;
vw=0.0;
for i=1:nx
    %Top
    f(i,ny,4)=f(i,ny,2);
    f(i,ny,7)=r*f(i,ny,5)+(1-r)*f(i,ny,6)+2*r*rho(i,ny)*w(7)*(cx(7)*uw+cy(7)*vw)/(2*c2);
    f(i,ny,8)=r*f(i,ny,6)+(1-r)*f(i,ny,5)+2*r*rho(i,ny)*w(8)*(cx(8)*uw+cy(8)*vw)/(2*c2);
    %Below
    f(i,1,2)=f(i,1,4);
    f(i,1,5)=r*f(i,1,7)+(1-r)*f(i,1,8);
    f(i,1,6)=r*f(i,1,8)+(1-r)*f(i,1,7);
end
%Left boundary, Periodic

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
% 
% for i=1:nx
%     for j=1:ny
%         for k=1:9
%             f2(k,i,j)=f(i,j,k);
%         end
%         m=mm*f2(:,i,j);
%         %计算平衡动量
%         m_eq(9,1)=rho(i,j);
%         m_eq(1,1)=rho(i,j)*(-2+3*(u(i,j)^2+v(i,j)^2));
%         m_eq(2,1)=rho(i,j)*(1-3*(u(i,j)^2+v(i,j)^2));
%         m_eq(3,1)=rho(i,j)*u(i,j);
%         m_eq(4,1)=-rho(i,j)*u(i,j);
%         m_eq(5,1)=rho(i,j)*v(i,j);
%         m_eq(6,1)=-rho(i,j)*v(i,j);
%         m_eq(7,1)=rho(i,j)*(u(i,j)^2-v(i,j)^2);
%         m_eq(8,1)=rho(i,j)*u(i,j)*v(i,j);
%         %施加外力，需处理u计算，guo书；
%         Fi(9,i,j)=0;
%         Fi(1,i,j)=6*(1-0.5*s1(1))*(u(i,j)*forcex+v(i,j)*forcey);
%         Fi(2,i,j)=-6*(1-0.5*s1(2))*(u(i,j)*forcex+v(i,j)*forcey);
%         Fi(3,i,j)=forcex;
%         Fi(4,i,j)=-(1-0.5*s1(4))*forcex;
%         Fi(5,i,j)=forcey;
%         Fi(6,i,j)=-(1-0.5*s1(6))*forcey;
%         Fi(7,i,j)=2*(1-0.5*s1(7))*(u(i,j)*forcex-v(i,j)*forcey);
%         Fi(8,i,j)=(1-0.5*s1(8))*(u(i,j)*forcey+v(i,j)*forcex);
%         %在动量空间松弛
%         m_temp=msl*(m-m_eq)-mminv*Fi(:,i,j);
%         for k=1:9
%             f2(k,i,j)=f2(k,i,j)-m_temp(k);
%             f(i,j,k)=f2(k,i,j);
%         end
%     end
% end
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