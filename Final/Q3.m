K=-1; 
M=1; 
N=1.5; 

Ubest = mean(U);
UStd=std(U); 

Vbest = mean(V);
VStd=std(V); 

Wbest = mean(W);
WStd=std(W); 

Q=(Ubest^K)*(Vbest^M)*(Wbest^N);
dQ=Q*sqrt((K*UStd/Ubest)^2+(M*VStd/Vbest)^2+(N*WStd/Wbest)^2);