
function  verticalwave(T)

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
B=zeros(xpt,ypt,ypt);
D=zeros(xpt,ypt,ypt);

p=1;
rho_1=0.01;
rho_2=0.010;

for i = 1:xpt
    for j = 1:ypt
        for k = 1:zpt
            X(i,j,k)=Lx*(i-1)/(xpt-1);
           Y(i,j,k)= Ly*(j-1)/(ypt);
           Z(i,j,k)= Lz*(k-1)/(zpt-1);
         
           
           
           U(i,j,k)= F(p,1);%+0.000000000000000000000000001*i;
           V(i,j,k)= F(p,2);%+0.000000000000000000000000001*i;
           W(i,j,k)= F(p,3);%+0.000000000000000000000000001*i;
           P(i,j,k)= F(p,4);%+0.000000000000000000000000001*i;
           B(i,j,k)= F(p,5);%+0.000000000000000000000000001*i;
           
           D(i,j,k)= 1-(0.001/1)*F(p,5)-tanh(Lz*((k-1)/(zpt-1)-0.5)/0.1);
           
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
p=reshape( P(:,yp,:), xpt,zpt);
b=reshape( B(:,yp,:), xpt,zpt);
d=reshape( D(:,yp,:), xpt,zpt);
Q=u.*u+v.*v+w.*w;



%contourf(x,z,p)
%contourf(x(:,1:end),z(:,1:end),u(:,1:end))
%quiver(x,z,u,w,'k')
%plot(w(1,:),z(1,:))


START1=1;
po=1;

subplot(2,1,1);
plot(u(po,:),z(po,:))
xlabel('U')
ylabel('z')

timename1='Vertical Velocity t= ';
timename2='s';
titletime = append(timename1, string(T*n*dt/100),timename2);
%title(titletime)
title('Vertical ')
axis([-0.4 0.4 0 Lz])

 set(gca,'YTick',0:Lz/2:Lz) 
AXLX=string(0:Lx/(pi*2):Lx/pi)+'pi';
 set(gca,'YTickLabel', string(AXLX))





subplot(2,1,2); 
plot(w(po,:),z(po,:))
xlabel('W')
ylabel('z')
grid off
axis([-0.15 0.15 0 Lz])

 set(gca,'YTick',0:Lz/2:Lz) 
AXLX=string(0:Lx/(pi*2):Lx/pi)+'pi';
 set(gca,'YTickLabel', string(AXLX))




end
