%{
[num,txt,all]=xlsread('Mohabir815214144.xls');
[~,c]=size(all);

for n=1:c;
ln=sprintf('%s=num(2:%d,%d);',txt{1,n},num(1,n)+1,n);
evalin('base',ln);
end
%}

%{

maxVal=max(A);
minVal=min(A);
meanVal=mean(A);
stdVal=std(A);

fprintf('max=%3.1f\n', maxVal);
fprintf('min=%3.1f\n', minVal);
fprintf('mean=%3.1f\n', meanVal);
fprintf('std=%3.1f\n', stdVal);

binsize=2;
[nn,bb]=hist(A(:),-32:binsize:68);

fprintf('area=%3.1f\n', sum(nn)*binsize);

normal_nn = nn/(sum(nn)*binsize);
fprintf('normalarea=%3.1f\n', sum(normal_nn)*binsize);

ideal_nn=(1/(stdVal*sqrt(2*pi))).*exp(-0.5*((bb-meanVal)./stdVal).^2);

hold on

bar(bb,normal_nn,1,'FaceColor','#80B3FF'); 
plot(bb,ideal_nn,'k','LineWidth', 1); 
xlabel('Value');                       % axis labels
ylabel('Normalized Count (P(A))');

legend('Normalized Data','Normal Distribution');

hold off
%}


maxValB=max(B);
minValB=min(B);
meanValB=mean(B);
stdValB=std(B);

fprintf('max=%3.1f\n', maxValB);
fprintf('min=%3.1f\n', minValB);
fprintf('mean=%3.1f\n', meanValB);
fprintf('std=%3.1f\n', stdValB);

[nnB,bbB]=hist(B(:),0:1:2);

fprintf('area=%3.1f\n', sum(nnB)*1);
normal_nnB = nnB/(sum(nnB)*1);

hold on

bar(bbB,normal_nnB,1,'FaceColor','#80B3FF'); 
plot(linspace(-0.5, 2.5, 10),(1/3)*ones(size(linspace(-0.5, 2.5, 10))),'k','LineWidth', 1); 

xlabel('Value');                       % axis labels
ylabel('Normalized Count (P(B))');

legend('Data', 'Uniform Distribution');

hold off

%{
maxValC=max(C);
minValC=min(C);
meanValC=mean(C);
stdValC=std(C);

fprintf('max=%3.1f\n', maxValC);
fprintf('min=%3.1f\n', minValC);
fprintf('mean=%3.1f\n', meanValC);
fprintf('std=%3.1f\n', stdValC);

[nnC,bbC]=hist(C(:),0:1:4);

fprintf('area=%3.1f\n', sum(nnC)*1);
normal_nnC = nnC/(sum(nnC)*1);
ideal_nnC=((meanValC).^bbC)*exp(-meanValC)./factorial(bbC);


hold on

bar(bbC,normal_nnC,1,'FaceColor','#80B3FF');  
plot(bbC,ideal_nnC,'k','LineWidth', 1); 

xlabel('Value');                       % axis labels
ylabel('Normalized Count (P(C))');

legend('Data', 'Poisson Distribution');

hold off
%}