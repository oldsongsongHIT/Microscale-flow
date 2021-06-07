%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear all; clc;
nx=101;ny=226;
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
%��ʼ��ù������
% [x,y] = meshgrid(1:nx,1:ny);%���ӿռ�
for i=1:nx
    for j=1:ny
        if (y1(j)>131.25+0.625*x1(i))
            obst(i,j)=2;
        elseif (y1(j)<31.25+0.625*x1(i))
            obst(i,j)=1;
        end
    end
end
obst(:,1)=0;
%��ù������
uo=0.00;  %�?般取0.1数量�?,�?般小�?0.3
alphag=1/3; %扩散系数
omegag=1.0./(3.*alphag+0.5);  %扩散系数对应松弛时间

count=0; cycle=40;%每隔cycle输出�?次结�?
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
% addpath(genpath('E:\LBM\Code\Reaction_flow'));
tic;
while error>tol
    [g]=gcol(nx,ny,u,v,cx,cy,omegag,g,rhog,w); %collision
    [g]=stream(g); %streaming
    [g]=gbound(nx,ny,w,g,rhog,alphag,omegag,obst); % boundary conditions
    %组分扩散-Obstacle
    % [g]=gobstc(nx,ny,g,opp,obst);
    rhog=sum(g,3);%%计算浓度
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
toc;
function [g]=gbound(nx,ny,w,g,rhog,alphag,omegag,obst)
%��ֶ������̵ı߽�����
%Boundary condition:
co=1.0;
% 
% for j=1:32
%     g(1,j,1)=g(2,j,1);
%     g(1,j,5)=g(2,j,5);
%     g(1,j,8)=g(2,j,8);
% end
% for j=1:94
%     g(nx,j,3)=g(nx-1,j,3);
%     g(nx,j,6)=g(nx-1,j,6);
%     g(nx,j,7)=g(nx-1,j,7);
% end
for j=33:132
g(1,j,1)=g(nx,j+62,1);
g(1,j,5)=g(nx,j+62,5);
g(1,j,8)=g(nx,j+62,8);
end
%right hand boundary,  zhouqi
for j=95:194
g(nx,j,3)=g(1,j-62,3);
g(nx,j,7)=g(1,j-62,7);
g(nx,j,6)=g(1,j-62,6);
end
%Top boundary, reaction
for i=1:nx
    for j=1:ny
        if obst(i,j)==1
            g(i,j,2)=-g(i,j,4);
            g(i,j,5)=-g(i,j,7);
            g(i,j,6)=-g(i,j,8);
        end
        if obst(i,j)==2
            g(i,j,4)=(w(2)+w(4))*co-g(i,j,2);
            g(i,j,8)=(w(6)+w(8))*co-g(i,j,6);
            g(i,j,7)=(w(5)+w(7))*co-g(i,j,5);
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
for i=2:nx-1
    for j=2:ny-1
        if obst(i,j)==1
            for k=1:9
                g(i,j,k) = g(i,j,opp(k));
            end
        elseif obst(i,j)==2
            % for k=9:-1:1
                % g(i,j,k) = g(i,j,opp(k));
            % end
			for k=1:9
                g(i,j,k) = g(i,j,opp(k));
            end
         end
    end
end
g(1,1,5)=g(1,1,7); g(nx,1,6)=g(nx,1,8); 
g(1,ny,8)=g(1,ny,6);g(nx,ny,7)=g(nx,ny,5);


end
function result(nx,ny,u,v,rhog,count,obst)
Cg=rhog;
pi=0;
xn=[1:nx]'; yn=[1:ny]';
for i=1:nx
    for j=1:ny
        if obst(i,j)==1||obst(i,j)==2
            Cg(i,j)=0.0;
        end
        pi=pi+1;
        peess(pi,:)=[xn(i),yn(j),u(i,j),v(i,j),Cg(i,j)];
    end
end
filename=['F:\LBM_code\date\' num2str(count) '-tecplot2d.dat'];
%   address是储存位置，这里的num2str是为了在循环输出dat数据文件中使用，如果只有�?个文件可以忽�?
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
