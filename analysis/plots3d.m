function  plots3d(T)

filename = './files/parameters.txt';
F = importdata(filename);

xpt=F(1,1);ypt=F(1,2);zpt=F(1,3);n  = F(1,4);dt =F(1,5);
Lx =F(2,1);Ly =F(2,2);Lz =F(2,3);Re = F(2,4);Fr =F(2,5);





na1='./files/U';
na2='.txt';
filename = append(na1, string(T),na2);
F = importdata(filename);

%xpt=ypt;zpt=ypt;
ypt=xpt;
           x= 0:1/(xpt-1):1;
           y= 0:1/(ypt-1):1;
           z= 0:1/(zpt-1):1;
         
U=zeros(xpt,ypt,ypt);
V=zeros(xpt,ypt,ypt);
W=zeros(xpt,ypt,ypt);
P=zeros(xpt,ypt,ypt);

p=1;

for i = 1:xpt
    for j = 1:ypt/2
        for k = 1:zpt

           
           
           U(i,j,k)= F(p,1);
           V(i,j,k)= F(p,2);
           W(i,j,k)= F(p,3);
           P(i,j,k)= F(p,4);
           B(i,j,k)= F(p,5);
           
           U(i,j+ypt/2,k)= F(p,1);
           V(i,j+ypt/2,k)= F(p,2);
           W(i,j+ypt/2,k)= F(p,3);
           P(i,j+ypt/2,k)= F(p,4);
           B(i,j+ypt/2,k)= F(p,5);
           
           p=p+1;
        end
    end
end

[X Y Z]=meshgrid(x,y,z);




yslice = [0  0.50  1];
%xslice = [ 0.25 0.50 0.75 ];

xslice = [];
zslice = [];
slice(X,Y,Z,B,xslice,yslice,zslice);
hold on
%quiver3(X,Y,Z,V,U,W,'k')
hold off
colormap(jet(64));
axis([0 1 0 1 0 1]);
maxu=max(max(max(U)));
minu=min(min(min(U)));
caxis([minu,maxu]);
colorbar('vertical');
shading interp;
colorbar
xlabel('y')
ylabel('x')
zlabel('z')



end

