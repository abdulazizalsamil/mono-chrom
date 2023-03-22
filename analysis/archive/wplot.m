function  wplot(T)


na1='./files/U';
na2='.txt';
filename = append(na1, string(T),na2);
F = importdata(filename);

xpt=F(1,1);ypt=F(1,2);zpt=F(1,3);n = F(1,4);
dt =F(2,1);Lx =F(2,2);Ly =F(2,3);Lz= F(2,4);
Ux =F(3,1);Vy =F(3,2);Wz =F(3,3);Ps= F(3,4);
Ts =F(4,1);rho=F(4,2);mu =F(4,3);Re = F(4,4);




F=F(5:end,:);


X=zeros(xpt,ypt,ypt);
Y=zeros(xpt,ypt,ypt);
Z=zeros(xpt,ypt,ypt);


U=zeros(xpt,ypt,ypt);
V=zeros(xpt,ypt,ypt);
W=zeros(xpt,ypt,ypt);
P=zeros(xpt,ypt,ypt);

p=1;
for i = 1:xpt
    for j = 1:ypt
        for k = 1:zpt
          
            X(i,j,k)= 1-(i-1)/(xpt)-1/xpt;
              if k>zpt/4
              % X(i,j,k)=1-(i-1)/(xpt)*((1)+1/4-(k-1)*(1/zpt))-1/xpt;
              end
           Y(i,j,k)= j/(ypt);
           Z(i,j,k)= k/(zpt);
         
           
           
           U(i,j,k)= F(p,1);
           V(i,j,k)= F(p,2);
           W(i,j,k)= F(p,3);
           P(i,j,k)= F(p,4);
           p=p+1;
        end
    end
end

yp=ypt/2;
if ypt==1
    yp=1;
    
end

x=reshape( X(:,yp,:), xpt,zpt);
z=reshape( Z(:,yp,:), xpt,zpt);

u=reshape( U(:,yp,:), xpt,zpt);
v=reshape( V(:,yp,:), xpt,zpt);
w=reshape( W(:,yp,:), xpt,zpt);

plot(w(xpt/2,:),z(xpt/2,:))
axis([-1 1 0 1])
xlabel("W")

ylabel("z")


end
