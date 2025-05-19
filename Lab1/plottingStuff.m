

x=1./((allL*0.001).^2);
y=allf;

x_s = [ones(size(x)) x];

w=(allL*0.001).^2;
W=diag(w);

m=x_s\y;

x_avg = mean(x);

y_avg = mean(y);

N = size(y,1);

r = sum((x-x_avg).*(y-y_avg))/sqrt(sum((x-x_avg).^2)*sum((y-y_avg).^2));

y_std = sqrt(sum((y-x_s*m).^2)./(N-2));

disp(y_std);

m_std = y_std*1/sqrt(N.*(sum(x.^2)./N - (sum(x)./N).^2));

disp(m);

disp(m_std);

hold on

scatter(x,y,'blue');
plot(x,x_s*m,'blue');
xlabel('1/L^2 (1/m^2)') 
ylabel('Freq (hz)') 

legend('Data','Fit')

hold off
