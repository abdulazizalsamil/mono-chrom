c=fix(clock);
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = './plots/M'+string(c(1))+string(c(2))+string(c(3))+string(c(4))+string(c(5))+string(c(6))+'.gif';
filen=100;


dt=0.001;
dy=0.01;
dz=0.01;

y = 0:dy:1;
z = 0:dz:1;
[Y,Z] = meshgrid(y,z);


k=1*2*pi;
n=10*2*pi;
Fr=0.01;
sigma=k/(Fr*(sqrt(k^2+n*2)));
N=0;
for T = 0:filen
N=N+10;
    v =(-n/k)*cos(k.*Y+n.*Z-sigma*dt*N);
    w = cos(k.*Y+n.*Z-sigma*dt*N);
  
    
    
tiledlayout(2,2)

nexttile
contourf(Y,Z,v)
title('U contour plot')
colorbar


nexttile
contourf(Y,Z,w)

title('W contour plot')
colorbar

nexttile
plot(Y(1,:),w(1,:))
hold on
plot(Y(1,:),v(1,:))
hold off
title('X W at z=0 plot')
xlabel("x")
ylabel("w")


nexttile
plot(w(:,50),Z(:,50))
hold on
plot(v(:,50),Z(:,50))
hold off
xlabel("W")
ylabel("z")
title('W Z at x=0.5 plot')

    
    
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if T == 0 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end
end