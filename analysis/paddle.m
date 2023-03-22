function  paddle(T)
zpt=64;
k=(1:zpt)-1;
dt=0.01;
u=exp((-(zpt-1)+k)*0.11)*cos(2*pi*dt*T);


plot(u,k/zpt)
axis([-1 1 0 1])
end

