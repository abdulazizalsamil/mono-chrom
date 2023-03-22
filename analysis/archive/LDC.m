T=100;
na1='./files/U';
na2='.txt';
filename = append(na1, string(T),na2);
F = importdata(filename);

xpt=F(1,1);ypt=F(1,2);zpt=F(1,3);n = F(1,4);
dt =F(2,1);Lx =F(2,2);Ly =F(2,3);Lz= F(2,4);
Ux =F(3,1);Vy =F(3,2);Wz =F(3,3);Ps= F(3,4);
Ts =F(4,1);rho=F(4,2);mu =F(4,3);Re = F(4,4);

G1 = importdata('ghiav.txt');
G=flipud(G1);
ghiax=G(:,1);
ghiav100=G(:,2);
ghiav400=G(:,3);
ghiav1000=G(:,4);

Gu1 = importdata('ghiau.txt');
Gu=flipud(Gu1);
ghiay=Gu(:,1);
ghiau100=Gu(:,2);
ghiau400=Gu(:,3);
ghiau1000=Gu(:,4);


F=F(3:end,:);


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
           X(i,j,k)= i/(xpt-1);
           Y(i,j,k)= j/(ypt-1);
           Z(i,j,k)= k/(zpt-1);
         
           
           
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

%u=-1*u;
%w=-1*w;
m=xpt/2;

tiledlayout(2,2)

nexttile
contourf(x,z,u)
xlabel('x')
ylabel('z')
colorbar
colormap(jet(64));
grid off
shading interp;


nexttile
contourf(x,z,w)
xlabel('x')
ylabel('z')
colorbar
colormap(jet(64));
grid off
shading interp;


nexttile
hold on
plot( u(m,:),(z(m,:)),'k');
plot( ghiau100,ghiay,'k*' );
grid on
xlabel('u')
ylabel('z')
axis([-0.25 1 0 1])

%legend('present results','Ghias results', 'Location','southeast','FontSize',15)
%legend('boxoff')

hold off


nexttile


hold on
plot( x(:,m),w(:,m),'k');
plot( ghiax,ghiav100,'k*' );
grid on
xlabel('x')
ylabel('w')
axis([0 1 -0.3 0.2])

%legend('present results','Ghias results', 'Location','northeast','FontSize',15)
%legend('boxoff')

hold off