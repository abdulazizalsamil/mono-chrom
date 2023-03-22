function  plots2d(T)


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
            X(i,j,k)=Lx*(i)/(xpt);
           Y(i,j,k)= (j)/(ypt);
           Z(i,j,k)= Lz*(k)/(zpt);
         
           
           
           U(i,j,k)= F(p,1);
           V(i,j,k)= F(p,2);
           W(i,j,k)= F(p,3);
           P(i,j,k)= F(p,4);
           B(i,j,k)= F(p,5);
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

subplot(3,1,1);
plot(x(:,zpt/2),u(:,zpt/2),".-")

 set(gca,'XTick',0:Lx/4:Lx) 
AXLX=string(0:Lx/(4*pi):Lx/(pi))+'pi';
 set(gca,'XTickLabel', string(AXLX))

timename1='time = ';
timename2=' s';
titletime = append(timename1, string(T*n*dt/100),timename2);
title(titletime)

axis([0 Lx -0.00001 0.00001])
xlabel('x')
ylabel('Horizontal Velocity')
grid on
%hold on



subplot(3,1,2);

plot(x(:,zpt/2),w(:,zpt/2),".-")
 set(gca,'XTick',0:Lx/4:Lx) 
AXLX=string(0:Lx/(4*pi):Lx/(pi))+'pi';
 set(gca,'XTickLabel', string(AXLX))

axis([0 Lx -0.0001 0.0001])
xlabel('x')
ylabel('Vertical Velocity')
grid on
%hold on


subplot(3,1,3);

plot(x(:,zpt/2),b(:,zpt/2),".-")

 set(gca,'XTick',0:Lx/4:Lx) 
AXLX=string(0:Lx/(4*pi):Lx/(pi))+'pi';
 set(gca,'XTickLabel', string(AXLX))
axis([0 Lx -0.0001 0.0002])
xlabel('x')
ylabel('Buoyancy')
grid on
%hold on
end
