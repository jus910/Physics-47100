function analysisfunc(filename)

M = readmatrix(filename);
t1=M(:,1);
v1=M(:,2);
t2=M(:,3);
v2=M(:,4); 



hold on

plot(t1,v1,'blue')
plot(t2,v2,'red')
xlabel('Time (s)') 
ylabel('Voltage (V)') 

legend('Piezo Force','Voltage due to Contact')

hold off