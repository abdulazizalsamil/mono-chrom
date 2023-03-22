T=3;


filename = './files/parameters.txt';
F = importdata(filename);

xpt=F(1,1);ypt=F(1,2);zpt=F(1,3);n  = F(1,4);dt =F(1,5);
Lx =F(2,1);Ly =F(2,2);Lz =F(2,3);Re = F(2,4);Fr =F(2,5);





na1='./files/U';
na2='.txt';
filename = append(na1, string(T),na2);
F = importdata(filename);



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
           Z(i,j,k)= 1*k/(zpt);
         
           
           
           U(i,j,k)= F(p,1);
           V(i,j,k)= F(p,2);
           W(i,j,k)= F(p,3);
           P(i,j,k)= F(p,4);
           p=p+1;
        end
    end
end

xp=xpt/2;
if xpt==1
    xp=1;
    
end

y=reshape( Y(xp,:,:), ypt,zpt);
z=reshape( Z(xp,:,:), ypt,zpt);

py=reshape( P(xp,:,:), ypt,zpt);
u=reshape( U(xp,:,:), ypt,zpt);
v=reshape( V(xp,:,:), ypt,zpt);
w=reshape( W(xp,:,:), ypt,zpt);
figure(1)
tiledlayout(2,2)

nexttile
contourf(y,z,v)
title('v contour plot at x midpoit')
colorbar


nexttile
contourf(y,z,w)
title('W contour plot at x midpoit')
colorbar

nexttile
plot(y(:,1),w(:,1))
title('Y W at z=0 plot')
xlabel("y")
ylabel("w")


nexttile
plot(w(ypt/2,:),z(ypt/2,:))
xlabel("W")
ylabel("z")
title('W Z at x=0.5 plot')
max(max(py))

figure(2)
tiledlayout(1,2)
nexttile
plot(w(ypt,:),z(ypt,:),'k')
hold on
plot(w(1,:),z(1,:))
hold off
title('W')


nexttile
plot(v(1,:),z(1,:),'k')
hold on
plot(v(ypt,:),z(ypt,:))
hold off
title('V')
