clear all
clc
xpt=16;
zpt=xpt;
dx=1/(xpt-1);
dz=1/(zpt-1);
X=zeros(xpt,zpt);
Z=zeros(xpt,zpt);
length=1/4;

A=1/(1-length);

for i=1:xpt
   for k=1:zpt
    
    X(i,k)=(i-1)*dx;
    Z(i,k)=(k-1)*dz;
    
    
    
    end
end
hold on
for i=1:xpt
plot(X(i,:),Z(i,:),'k')
  
end

for k=1:zpt
plot(X(:,k),Z(:,k),'k')  
end


hold on



%This is for velocity
for i=0:xpt-1

   for k=0:zpt-1-ceil(A*i)
    
    plot(X(i+1,k+1),Z(i+1,k+1),'ko')
   end
   
end

%pressure
for i=0:xpt-1%-fix((xpt+0.1)*length)

   k=zpt-1-ceil(A*(i+0.1));
   if k>-0.1 && length>0
   plot(X(i+1,k+1),Z(i+1,k+1),'r*')
   end
   
  
   
      k=(zpt-1-ceil(A*(i)));
      
      
   if k>-0.1   && length>0
   plot(X(i+1,k+1),Z(i+1,k+1),'r*')
   end
   
   
   
end



%for the derivative
  for i=1:xpt-1
      
      
   for k=0:zpt-1
     if k>=0 && k==zpt-1+1-ceil(A*(i))
     plot(X(i+1,k+1),Z(i+1,k+1),'b^')
     end  
       
   end
      
      
      
             k=zpt-ceil(A*(i));
   if k>=0
   %plot(X(i+1,k+1),Z(i+1,k+1),'b^')
   end  
      
  end
   

s=128;
S=0:s;
X1=S./s;
z=1-A*(X1);
axis([0 1 0 1])
plot(X1,z,'k') 
hold off
