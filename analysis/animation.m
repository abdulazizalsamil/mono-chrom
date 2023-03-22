
c=fix(clock);
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = './plots/M'+string(c(1))+string(c(2))+string(c(3))+string(c(4))+string(c(5))+string(c(6))+'.gif';
filbe=0;
filen=200;

global minU maxU minV maxV minW maxW minP maxP minB maxB minD maxD  minE maxE

           minU=0;maxU=0;minV=0;maxV=0;minW=0;maxW=0;
           minP=0;maxP=0;minB=0;maxB=0;minD=1;maxD=1;
           minE=1;maxE=1;


[minU,maxU,minV,maxV,minW,maxW,minP,maxP,minB,maxB,minD,maxD,minE,maxE]=MaxValue(filbe,filen);

for T = filbe:1:filen
    %yplot4(T)
    %plot4(T)
    %plots2d(T)
    %plots3d(T)
    %wplot(T)
    %paddle(T)
    %wb_bc(T)
    wavemaker(T)
    %verticalwave(T)
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if T == filbe 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end
end

