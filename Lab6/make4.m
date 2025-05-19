
binsize = 0.01; 

[nny,bby]=hist(vys(:),-.5:binsize:.5);
[nnx,bbx]=hist(vxs(:),-.5:binsize:.5);

meanX = mean(vxs(:));
stdX = std(vxs(:));

meanY = mean(vys(:));
stdY = std(vys(:));

ideal_nnx=(1/(stdX*sqrt(2*pi))).*exp(-0.5*((bbx-meanX)./stdX).^2);
normal_nnx = nnx/(sum(nnx)*binsize);
normal_nny = nny/(sum(nny)*binsize);


hold on

bar(bbx,normal_nnx,1,'red'); 
bar(bby,normal_nny,1,'FaceColor','#80B3FF'); 
plot(bbx,ideal_nnx,'k','LineWidth', 1.5); 
xlabel('Velocity');                       % axis labels
ylabel('Normalized Count');

legend('Normalized Vx','Normalized Vy', 'Ideal Normal Distribution Vx' )

hold off
disp(meanX);
disp(stdX);
disp(sum(normal_nnx)*binsize);