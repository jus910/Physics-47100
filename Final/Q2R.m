

x=log(Rx);
y=log(Ry);

x_s = [ones(size(x)) x];
m=x_s\y;

x_avg = mean(x);

y_avg = mean(y);

N = size(y,1);

y_std = sqrt(sum((y-x_s*m).^2)./(N-2));

r = sum((x-x_avg).*(y-y_avg))/sqrt(sum((x-x_avg).^2)*sum((y-y_avg).^2));

m_std = y_std*1/sqrt(N*(sum(x.^2)/N - (sum(x)/N)^2));

b_std = y_std*sqrt((sum(x.^2)./N))/sqrt(N.*(sum(x.^2)./N - (sum(x)./N).^2));


hold on

scatter(x,y); 
plot(x,x_s*m,'k'); 
plot(x,x_s*m+y_std,'--k'); 
plot(x,x_s*m-y_std,'--k'); 
xlabel('log(Rx)');                       % axis labels
ylabel('log(Ry)');

legend('Data','Fit')

hold off
