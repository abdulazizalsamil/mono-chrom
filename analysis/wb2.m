clear all
clc
T=0;

zpt=128;xpt=128;


dt=0.05;

X=zeros(xpt,zpt);Z=zeros(xpt,zpt);
W=zeros(xpt,zpt);U=zeros(xpt,zpt);
Lx=100;Lz=1;
dx=Lx/(xpt);dz=Lz/(zpt);


%const
N1=32;
x_b=Lx/2;z_0=Lz/2;
lampda_bar=0.9*(Lz/2);
epsinon=0.01;
g=9.8;
sigma_bar=epsinon*g*(2*pi)/(lampda_bar);
delta_sigma=(2/3)*(sigma_bar/(N1-1));
k_bar=2*pi/lampda_bar;
cg=0.5*(sigma_bar)/(k_bar);
H=1;
for i = 1:xpt
    for k = 1:zpt
           X(i,k)= i*dx;
           Z(i,k)= k*dz;
          
           
     end
end
sigma_j=0;
k_j=0;
for N = 1:N1
    for k = 1:zpt
        sigma_j=k*sigma_bar; k_j=k_bar;
       W(1,k)=0.01*log(cosh(dz*k));%W(1,k)+exp(-((k-1)*dz/H)^2)*( sigma_j*( sin( (sigma_j/cg-k_j)*x_b-sigma_j*dt*T ) )  );
           
     end
end



%quiver(X,Z,U,W,'k')
plot(W(1,:),Z(1,:))
%axis([(-0.01*Lx+min(min(X))) max(max(X)) min(min(Z)) max(max(Z)) ])
xlabel('non-constant part of b')
ylabel('z')