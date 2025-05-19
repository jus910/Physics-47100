tic;
override = {'KE0', 1, 'phi_set', 0.7};
mdNVE000;
toc;
clear('override');

t=(1:Nt)*dt;                         % make a time variable from dt to TT
plot(t,[Ek Ep Ek+Ep],'linewidth',3)  % concatenate variable to add to plot
legend('Kinetic Energy','Potential Energy','Total Energy');
xlabel('Time')                       % axis labels
ylabel('Energy')