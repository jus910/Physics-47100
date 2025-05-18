function testingFit()

filename = "testing6.csv";

M = readmatrix(filename);
y=M(:,2);
x=M(:,1);
w=M(:,3);
x_one=ones(size(x));

x_s = [x_one x];
disp(x_s);

const = x_s \y;

W=diag(w);
constW=x_s'*W*x_s\(x_s'*W*y);
disp(constW);

hold on

plot(x,y,'o','MarkerSize',3,'Color','#7E2F8E')
plot(x,x_s*constW,'-','MarkerSize',3,'Color','red')
plot(x,x_s*const,'-','MarkerSize',3,'Color','blue')
%xline(17,'red')
xlabel('1/r^2 (1/m^2)') 
ylabel('Counts') 

%legend('Data','Time of Restriction')

hold off