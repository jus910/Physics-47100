%%  Setup
N=1e4;            % Number of points
t=0:.01:2*pi;     % helper list

rng(sum('Your Name'),'twister'); % reset random number generator

%% Monte Carlo integration for the area of a circle

u=2*rand(2,N)-1;  % Choose random point in the 2D plane [-1 1]x[-1 1]
x=u(1,:);         % Nicer names
y=u(2,:);

ii=(x.^2+y.^2<=1); % True for points inside circle of radius 1

%% Plot points
plot(x(ii),y(ii),'b.',x(~ii),y(~ii),'r.',sin(t),cos(t),'k'); 
axis('square')
set(gca,'fontsize',20);
xlabel('x-axis');
ylabel('y-axis');
title('Random Points');
print -depsc2 random.eps

%% Calculate and display Pi

myPI=4*sum(ii)/N;
disp(myPI)
disp(sum(ii))