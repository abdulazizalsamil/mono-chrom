function  wb_bc(T)



zpt=128;xpt=128;


dt=0.01;
Et=10000*dt

X=zeros(xpt,zpt);Z=zeros(xpt,zpt);
W=zeros(xpt,zpt);U=zeros(xpt,zpt);
Lx=100;Lz=1;
dx=Lx/(xpt);dz=Lz/(zpt);


%const
trwidth=0.1;
Fr=1.0;
N1=32;
x_b=Lx/2;z_0=Lz/2;
lampda_bar=0.9*(2*pi/2);
epsinon=0.01;
g=1.0;

sigma_bar=Fr^2*trwidth*(2*pi)/(lampda_bar);

delta_sigma=(2/3)*(sigma_bar/(N1-1));
k_bar=2*pi/lampda_bar;
cg=0.5*(sigma_bar)/(k_bar);
H=0.2;
sigma_j=lampda_bar/20;
k_j=sigma_j.^2/(Fr^2*trwidth);

Bz = 0;
DBz= 0;

t1=2*(k_j/sigma_j)*x_b;

for i = 1:xpt
    for k = 1:zpt
           X(i,k)= i*dx;
           Z(i,k)= k*dz;
          
           
     end
end
for N = 1:N1
        
        
    for k = 2:zpt-1
        Bz = exp(-  ( ((k-1)*dz-0.5)/H)^2  );
        %DBz=-(((k-1)*dz-0.5)*(2.0/(H*H)))*exp(-  ( ((k-1)*dz-0.5)/H)^2  );
        DBz=-(((k-1)*dz-0.5)*(2.0/(H*H*1)))*exp(-  ( ((k-1)*dz-0.5)/H)^2  );
           
       W(1,k)=W(1,k)+Bz*( sigma_j*( sin( (sigma_j/cg-k_j)*x_b-sigma_j*dt*T ) )  );
       U(1,k)=U(1,k)+DBz*((sigma_j/k_j)*( cos( (sigma_j/cg-k_j)*x_b-sigma_j*dt*T ) )  );
          
    end
        sigma_j=sigma_j+delta_sigma;
        k_j=sigma_j.^2/(Fr^2*trwidth);
end

t2=2*(k_j/sigma_j)*x_b-t1

%quiver(X,Z,U,W,'k')

%figure(1)
plot(W(1,:),Z(1,:),'b')
legend('W')
grid on

xlabel('Velocity')
ylabel('Z')

%title1=compose("Lx=%1.1f,Lz=%1.1f,H=%1.1f",Lx,Lz,H);
title1=compose("W");
title(title1)


figure(2)
plot(U(1,:),Z(1,:),'r')
legend('U')
grid on
xlabel('Velocity')
ylabel('Z')
title('U')

%t_2=(k_j/sigma_j)*x_b
%k_j
%sigma_j
%title('B(z)=sin(2*pi*Z)')
%axis([(-0.01*Lx+min(min(X))) max(max(X)) min(min(Z)) max(max(Z)) ])
%axis([-100 100 0 Lz])
end
