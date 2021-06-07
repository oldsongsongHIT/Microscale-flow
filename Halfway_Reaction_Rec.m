%LBM- 2-D2Q9, Heated channel, note that c2=1/3, w9=4/9,
% w1-4=1/9, and w5-w8, 1/36
clear all; 
nx=100;ny=80;
Maxcount=10;
gpre=zeros(nx,ny,9);g=zeros(nx,ny,9);rhog=ones(nx,ny);
u=zeros(nx,ny);v=zeros(nx,ny);obst=zeros(nx,ny);

w=[1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36 4/9];
cx = [1 0 -1 0 1 -1 -1 1 0];
cy = [0 1 0 -1 1 1 -1 -1 0];
opp = [3, 4, 1, 2, 7, 8, 5, 6, 9];
c2=1./3.;
dx=1.0;dy=1.0;
x1=(0:1:nx); y1=(0:1:ny);
alphag=1/6; %鎵╂暎绯绘暟
omegag=1.0./(3.*alphag+0.5);  %鎵╂暎绯绘暟瀵瑰簲鏉惧紱鏃堕棿
rhog(:)=0.0;
%rection
count=0; cycle=40;%姣忛殧cycle杈撳嚭涓?娆＄粨鏋?
tol=1.0e-8; error=10.;erso=0.0;
%Main Loop
% addpath(genpath('E:\LBM\Code\Reaction_flow'));
tic;
while error>tol
    [g]=gcol(nx,ny,u,v,cx,cy,omegag,g,rhog,w); %collision
    gpre=g;
    [g]=stream(g); %streaming
    [g]=gbound(nx,ny,w,g,gpre,rhog,alphag,omegag,opp); % boundary conditions    
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
%         result(nx,ny,u,v,rhog,count,obst);
%     end
    count=count+1;
    if count>Maxcount
        break
    end
end
result(nx,ny,u,v,rhog,count,obst);
toc;
function [g]=gbound(nx,ny,w,g,gpre,rhog,alphag,omegag,opp)
%组分对流方程的边界条件
%Boundary condition:
co=10.0; %left hand all temperature
kr=0.1;
%left boundary, inlet temp. co
for j=1:ny
g(1,j,1)=(w(1)+w(3))*co-gpre(1,j,3);
g(1,j,5)=(w(5)+w(7))*co-gpre(1,j,7);
g(1,j,8)=(w(6)+w(8))*co-gpre(1,j,6);
end
% 自由出流边界条件
g(nx,:,3)=g(nx-1,:,3);
g(nx,:,7)=g(nx-1,:,7);
g(nx,:,6)=g(nx-1,:,6);
%bottom boundary, ad:abatic

%Top boundary, reaction
ceq=1.0;
for i=1:nx
    j=ny;
    b1=alphag;  b2=kr;  b3=kr*ceq;
    cw=(rhog(i,j)+0.5*b3/b1)/(1.0+0.5*b2/b1);
    g(i,j,4)=(w(2)+w(4))*cw-gpre(i,j,2);
    g(i,j,7)=(w(5)+w(7))*cw-gpre(i,j,5);
    g(i,j,8)=(w(6)+w(8))*cw-gpre(i,j,6);
end
g(:,1,2)=g(:,2,2);
g(:,1,5)=g(:,2,5);
g(:,1,6)=g(:,2,6);
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
function result(nx,ny,u,v,rhog,count,obst)
[X,Y]=meshgrid(0:nx, 0:ny);
filemat=strcat(num2str(count),'.mat');
save(filemat,'rhog');
for i=1:nx
    for j=1:ny
      if obst(i,j)==1
          rhog(i,j)=0;
      end
    end
end
filename=[num2str(count) '-2d.dat'];
%   address是储存位置，这里的num2str是为了在循环输出dat数据文件中使用，如果只有一个文件可以忽略
fid=fopen(filename,'wt');
fprintf(fid, 'VARIABLES = "X", "Y", "C"\n');
fprintf(fid, 'ZONE I=%d, J=%d, DATAPACKING=BLOCK, VARLOCATION=([3]=CELLCENTERED)\n', nx+1,ny+1);
fprintf(fid,'SOLUTIONTIME=%d\n',count);
for i=1:ny+1
fprintf(fid, '%g\t', X(i,:));
fprintf(fid, '\n');
end
fprintf(fid, '\n');
for j=1:ny+1
fprintf(fid, '%g\t', Y(j,:));
fprintf(fid, '\n');
end
fprintf(fid, '\n');
for i=1:ny
fprintf(fid, '%g\t', rhog(:,i));
fprintf(fid, '\n');
end
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
