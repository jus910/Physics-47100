function analysisfunc(filename)

M = readmatrix(filename);
x=M(:,1);
y=M(:,2);



hold on

plot(x,y,'-o','MarkerSize',5,'Color','#0072BD')
%yline(-6.28318530718,'red')
xlabel('x') 
ylabel('y') 

hold off