%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear all; clc;
nx=101;ny=81;
Maxcount=400000;
f=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=ones(nx,ny);
u=zeros(nx,ny);v=zeros(nx,ny);obst=zeros(nx,ny);
rho=ones(nx,ny);Cg=zeros(nx,ny);
w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];
c2=1./3.;
dx=1.0;dy=1.0;
x1=(0:1:nx); y1=(0:1:ny);
alphag=1/6; %鎵╂暎绯绘暟
omegag=1.0./(3.*alphag+0.5);  %鎵╂暎绯绘暟瀵瑰簲鏉惧紱鏃堕棿

%rection
count=0; cycle=40;%姣忛殧cycle杈撳嚭涓?娆＄粨鏋?
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
% addpath(genpath('E:\LBM\Code\Reaction_flow'));
tic;
while error>tol
    [g]=gcol(nx,ny,u,v,cx,cy,omegag,g,rhog,w); %collision
    [g]=stream(g); %streaming
    [g]=gbound(nx,ny,w,g,rhog,alphag,omegag,opp); % boundary conditions    
    rhog=sum(g,3);%%璁＄畻娴撳害
    ers=0.;
    for i =1:nx
        for j=1:ny
            ers=ers+u(i,j)*u(i,j)+v(i,j)*v(i,j)+rhog(i,j)*rhog(i,j);
        end
    end
    error=abs(ers-erso)/ers
    erso=ers;
    %杈撳嚭鏂囦欢
%     if mod(count,cycle)==0
%         result(nx,ny,u,v,rhog,count,obst,Cg);
%     end
    count=count+1;
    if count>Maxcount
        break
    end
end
result(nx,ny,u,v,rhog,count,obst,Cg);
toc;
function [g]=gbound(nx,ny,w,g,rhog,alphag,omegag,opp)
%组分对流方程的边界条件
%Boundary condition:
co=10.0; %left hand all temperature
kr=0.1;
%left boundary, inlet temp. co
g(1,:,1)=(w(1)+w(3))*co-g(1,:,3);
g(1,:,5)=(w(5)+w(7))*co-g(1,:,7);
g(1,:,8)=(w(6)+w(8))*co-g(1,:,6);
%回弹边界，不控制角点
% for j=1:ny
% g(nx,j,3)=g(nx,j,1);
% g(nx,j,7)=g(nx,j,5);
% g(nx,j,6)=g(nx,j,8);
% end
% %bottom boundary, adiabatic
% for i=1:nx-1
% g(i,1,2)=g(i,1,4);
% g(i,1,5)=g(i,1,7);
% g(i,1,6)=g(i,1,8);

% % 回弹边界，控制角点
% for j=2:ny-1
% g(nx,j,3)=g(nx,j,1);
% g(nx,j,7)=g(nx,j,5);
% g(nx,j,6)=g(nx,j,8);
% end
% %右下角
% g(nx,1,6)=g(nx,1,8);
% g(nx,1,2)=g(nx,1,4);
% g(nx,1,3)=g(nx,1,1);
% % g(nx,1,5)=g(nx,1,7);
% % g(nx,1,7)=g(nx,1,5);
% %右上角
% g(nx,ny,3)=g(nx,ny,1);
% g(nx,ny,4)=g(nx,ny,2);
% % g(nx,ny,6)=g(nx,ny,8);
% g(nx,ny,7)=g(nx,ny,5);
% % g(nx,ny,8)=g(nx,ny,6);
% 
% %bottom boundary, adiabatic
% for i=1:nx-1
% g(i,1,2)=g(i,1,4);
% g(i,1,5)=g(i,1,7);
% g(i,1,6)=g(i,1,8);
% end
% 自由出流边界条件
g(nx,:,3)=g(nx-1,:,3);
g(nx,:,7)=g(nx-1,:,7);
g(nx,:,6)=g(nx-1,:,6);
%bottom boundary, ad:abatic
g(:,1,2)=g(:,2,2);
g(:,1,5)=g(:,2,5);
g(:,1,6)=g(:,2,6);
%Top boundary, reaction
ceq=1.0;
for i=1:nx
       C(i,ny)=(alphag.*rhog(i,ny-1)+kr*ceq)/(kr+alphag);      
%      g(i,ny,4)=((w(2)+w(4))*rhog(i,ny)-kr*(rhog(i,ny)-ceq))./2;
%      g(i,ny,7)=((w(5)+w(7))*rhog(i,ny)-kr*(rhog(i,ny)-ceq))./2;
%      g(i,ny,8)=((w(6)+w(8))*rhog(i,ny)-kr*(rhog(i,ny)-ceq))./2;
    g(i,ny,4)=(w(2)+w(4))*C(i,ny)-g(i,ny,2);
    g(i,ny,7)=(w(5)+w(7))*C(i,ny)-g(i,ny,5);
    g(i,ny,8)=(w(6)+w(8))*C(i,ny)-g(i,ny,6); 
end
end
% 缁勫垎瀵规祦鏂圭▼鐨凜ollition
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
% Obsticale-缁勫垎鎵╂暎鏂圭▼闅滅
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
function result(nx,ny,u,v,rhog,count,obst,Cg)
Cg=rhog;
pi=0;
xn=[1:nx]'-1; yn=[1:ny]'-1;
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            Cg(i,j)=0.0;
        end
        pi=pi+1;
        peess(pi,:)=[xn(i),yn(j),u(i,j),v(i,j),Cg(i,j)];
    end
end
filename=['F:\date\' num2str(count) '-tecplot2d.dat'];
fid=fopen(filename,'wt');
fprintf(fid,'variables= "x", "y", "U", "V", "C"\r\n');
fprintf(fid,'zone t="Frame 0"i=%d,j=%d,f=point\r\n',ny,nx);
fprintf(fid,'SOLUTIONTIME=%d\r\n',count);
fprintf(fid,'%8.4f %8.4f %8.4f %8.4f %8.4f \r\n',peess');
fclose(fid);
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
