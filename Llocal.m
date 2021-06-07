function [Loc] = Llocal(obst,nx,ny,Loc)
%计算局部特征参数；
%1.计算各流体点的最大园；
%2. 删除多余的圆
Rad=zeros(nx,ny);
%以内部点为圆心，Rad对应最大圆的半径
for i=2:nx-1
    for j=2:ny-1
        if obst(i,j)==0
            %该点与四边界距离
            dis(i,j,1)=nx*ny;%(i-0)^2;
            dis(i,j,2)=nx*ny;%(i-nx-1)^2;
            dis(i,j,3)=(j-0)^2;
            dis(i,j,4)=(j-ny-1)^2;
            %检测该点与固体障碍物点的距离
            for m=1:nx
                for n=1:ny
                    if obst(m,n)==1
                        k=5;
                        dis(i,j,k)=(i-m)^2+(j-n)^2;
                        k=k+1;
                    end 
                end
            end
            Rad(i,j)=min(dis(i,j,:));
        end
    end
end
Cir=Rad.^0.5;
%s删除多余圆,Cir
for i=1:nx
    for j=1:ny
        if Cir(i,j)~=0
            for m=1:nx
                for n=1:ny
                    if Cir(m,n)~=0&&(m~=i||n~=j)&&(i-m)^2+(j-n)^2<=(Cir(i,j)-Cir(m,n))^2
                        if Cir(i,j)>=Cir(m,n)
                        Cir(m,n)=0;
                        else
                            Cir(i,j)=0;
                        end
                    end
                end
            end
        end
    end
end
dis2=(nx*ny)^2*ones(nx,ny,nx,ny);
%求距离最近的内切圆的直径为特征长度
for i=1:nx
    for j=1:ny
        for m=1:nx
            for n=1:ny
                if Cir(m,n)~=0
                    dis2(i,j,m,n)=(i-m)^2+(j-n)^2;
                end
            end
        end
        dis3=squeeze(dis2(i,j,:,:));
        e=min(min(dis3));
        [row,vol]=find(dis3==e,1);
        Loc(i,j)=Cir(row,vol)*2;
    end
end
end
        
            

