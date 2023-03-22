
function [minU,maxU,minV,maxV,minW,maxW,minP,maxP,minB,maxB,minD,maxD,minE,maxE]=MaxValue(filbe,filen)

           minU=0;maxU=0;minV=0;maxV=0;minW=0;maxW=0;
           minP=0;maxP=0;minB=0;maxB=0;minD=0;maxD=0;
           
           minE=0;maxE=0;

           
filename = './files/parameters.txt';
F = importdata(filename);

xpt=F(1,1);ypt=F(1,2);zpt=F(1,3);n  = F(1,4);dt =F(1,5);
Lx =F(2,1);Ly =F(2,2);Lz =F(2,3);Re = F(2,4);Fr =F(2,5);

      
           
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
na1='./files/U';
na2='.txt';
for T = filbe:filen
    






filename = append(na1, string(T),na2);
F = importdata(filename);



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
           
           D(i,j,k)= -(0.01/0.1)*B(i,j,k)+1-0.01*tanh(Lz*((k-1)/(zpt-1)-0.5)/0.1);
           p=p+1;
        end
    end
end

if maxU<max(max(max(U)))
 maxU=max(max(max(U)));   
end 
if minU>min(min(min(U)))
 minU=min(min(min(U)));   
end 


if maxV<max(max(max(V)))
 maxU=max(max(max(U)));   
end 
if minV>min(min(min(V)))
 minV=min(min(min(V)));   
end 


if maxW<max(max(max(W)))
 maxW=max(max(max(W)));   
end 
if minW>min(min(min(W)))
 minW=min(min(min(W)));   
end 



if maxP<max(max(max(P)))
 maxP=max(max(max(P)));   
end 
if minP>min(min(min(P)))
 minP=min(min(min(P)));   
end 


if maxB<max(max(max(B)))
 maxB=max(max(max(B)));   
end 
if minB>min(min(min(B)))
 minB=min(min(min(B)));   
end 


if maxD<max(max(max(D)))
 maxD=max(max(max(D)));   
end 
if minD>min(min(min(D)))
 minD=min(min(min(D)));   
end 

if maxE<max(max(max(E)))
 maxE=max(max(max(E)));   
end 
if minE>min(min(min(E)))
 minE=min(min(min(E)));   
end 


end

end

