function  wavemaker(T)
global minU maxU   minW maxW   minB maxB  minE maxE

filename = './files/parameters.txt';
F = importdata(filename);

xpt=F(1,1);ypt=F(1,2);zpt=F(1,3);n  = F(1,4);dt =F(1,5);
Lx =F(2,1);Ly =F(2,2);Lz =F(2,3);Re = F(2,4);Fr =F(2,5);
LN =F(3,1);



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

E=zeros(xpt,ypt,ypt);

p=1;

if T==0
   jjjj=1; 
else
   jjjj=0;  
end

for i = 1:xpt
    for j = 1:ypt
        for k = 1:zpt
            X(i,j,k)=Lx*(i-1)/(xpt-1);
           Y(i,j,k)= Ly*(j-1)/(ypt);
           Z(i,j,k)= Lz*(k-1)/(zpt-1);
         
           
           
           U(i,j,k)= F(p,1)+0.000000000000000000000000001*i*jjjj;
           V(i,j,k)= F(p,2)+0.000000000000000000000000001*i*jjjj;
           W(i,j,k)= F(p,3)+0.000000000000000000000000001*i*jjjj;
           P(i,j,k)= F(p,4)+0.000000000000000000000000001*i*jjjj;
           B(i,j,k)= F(p,5)+0.000000000000000000000000001*i*jjjj;
           
           E(i,j,k)= U(i,j,k)*W(i,j,k);
           
           D(i,j,k)= (1.0/LN)*B(i,j,k)-tanh(Lz*((k-1)/(zpt-1)-0.5)/LN);
           p=p+1;
        end
    end
end


if T==0
   U(1,1,1)=1;
   W(1,1,1)=1;
  
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


subplot(4,1,1);
contourf(x(START1:end,:),z(START1:end,:),Q(START1:end,:),'edgecolor','none')
%hold on
%quiver(x(START1:end,:),z(START1:end,:),u(START1:end,:),w(START1:end,:),'k')
%hold off
xlabel('x')
ylabel('z')
colorbar
colormap(jet(64));
grid off
shading interp;
axis([0 Lx 0 Lz])
timename1='Energy at t= ';
timename2=' ';
titletime = append(timename1, string(T*n*dt/100),timename2);
title(titletime)
energymax =maxE;%max([maxW*maxW,maxU*maxU,maxU*maxW,maxU,maxW]);
energymin =minE;%min([minW*minW,minU*minU,minU*minW,minU,minW]);
caxis([energymin energymax])





subplot(4,1,2); 
contourf(x(START1:end,:),z(START1:end,:),u(START1:end,:),'edgecolor','none')
xlabel('x')
ylabel('z')
colorbar
colormap(jet(64));
grid off
shading interp;
axis([0 Lx 0 Lz])
title('horizontal velocity')
caxis([minU maxU])





subplot(4,1,3); 
contourf(x(START1:end,:),z(START1:end,:),w(START1:end,:),'edgecolor','none')
xlabel('x')
ylabel('z')
colorbar
colormap(jet(64));
grid off
shading interp;
axis([0 Lx 0 Lz])
title('vertical velocity')
caxis([minW maxW])

subplot(4,1,4); 
contourf(x(START1:end,:),z(START1:end,:),b(START1:end,:),'edgecolor','none')
xlabel('x')
ylabel('z')
colorbar
colormap(jet(64));
grid off
shading interp;
axis([0 Lx 0 Lz])
title('buoyance')
caxis([minB maxB])

end