


Qx_s = [ones(size(Qx)) Qx];
mQ=Qx_s\Qy;

Qx_avg = mean(Qx);

Qy_avg = mean(Qy);

N = size(Qy,1);

Qy_std = sqrt(sum((Qy-Qx_s*mQ).^2)./(N-2));

Qr = sum((Qx-Qx_avg).*(Qy-Qy_avg))/sqrt(sum((Qx-Qx_avg).^2)*sum((Qy-Qy_avg).^2));

Qm_std = Qy_std*1/sqrt(N*(sum(Qx.^2)/N - (sum(Qx)/N)^2));

Qb_std = Qy_std*sqrt((sum(Qx.^2)./N))/sqrt(N.*(sum(Qx.^2)./N - (sum(Qx)./N).^2));


hold on

scatter(Qx,Qy); 
plot(Qx,Qx_s*mQ,'k'); 
plot(Qx,Qx_s*mQ+Qy_std,'--k'); 
plot(Qx,Qx_s*mQ-Qy_std,'--k'); 
xlabel('Qx');                       % axis labels
ylabel('Qy');

legend('Data','Fit')

hold off