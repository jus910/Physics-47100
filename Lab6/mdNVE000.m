%% Molecular Dynamics Simulator
% <md.m> Mark D. Shattuck 7/22/2010

% revision history:
% 7/22/2010 Mark D. Shattuck <mds> md.m
%           MD Demo for HandsOn 2010
%           000 mds Initial conditions and visualization
% 7/24/2010 001 mds Add Euler
% 7/24/2010 002 mds Add Interaction detection and Force Law
% 7/24/2010 003 mds Add KE, PE, Nplotskip, and plotit
% 1/21/2011 010 mds Add Velocity Verlet
% 3/28/2020 000 mds convert NPE simulation for phys 471 mdNPE.m
% 4/18/2020 000 mds clean up add comments fix override
% 4/20/2020 000 mds convert NVE simulation for phys 471 mdNVE.m


%% Experimental Parmeters (overrideable see below)
N=80;           % Number of particles 
Dn=2;           % diameter of particles [can be 1xN list]
K=100;          % spring constant particle-particle
M=3;            % mass of particles [can be 1xN list]
B=0;            % linear velocity particle drag

KE0=3;          % Initial particles Kinetic energy

% initial condition type choose one
ic='hexagonal';               % hexagonal template
% ic='square';                  % square template
% ic='load';                    % load from file
% icfile.name='compcry1';       % load ic file name
% icfile.var={'x','y','ytp'}    % variables to load

% boundaries 
setDensity=true;  % density determines box size
phi_set=.5;       % used to set density phi=pi*sum(Dn.^2)/4/Lx/Ly

Gb=1;             % aspect ratio of box Gb=Lx/Ly

L0_set=20;        % if setDensity is false use Ly=max(Dn)*L0; Lx=Gb*Ly;        
   
TT=150;           % total simulation time

saveall=true;  % save all positions and velocities?
clean=true;    % final processing?
pauseit=false; % pause after initial conditions

%% Display Parameters
plotit=true;     % plot?
Nplotskip=150;   % number of timesteps to skip before plotting

%% Overrides
% when running this script any variable above this section can be overriden
% by defining a varible override={'var_to_change',new_value}; For example
% to set the energy E to 2.3 and initial condition ic='square' and
% plotit=false use:
%
% override={'KE0',2.3,'ic','square','plotit',false};
% mdNVE000;
%
% To use the defaults: clear('override');

if(exist('override','var'))   % does override exist?
  Nover=length(override);     
  if(mod(Nover,2)~=0)         % must be pairs
    error('Error: override must contain pairs!');
  else
    Nover=Nover/2;            % number of pairs
    for n=1:Nover;
      m=2*n-1;
      assignin('base',override{m},override{m+1}); % set values
    end
  end
end

%% Simulation Parmeters
dt=1e-2;           % integration time step
Nt=fix(TT/dt);     % number of time steps

%% Calculated values
if(numel(Dn)==1);  % convert to list if not already
  Dn=Dn*ones(1,N); 
end
if(numel(M)==1);
  M=M*ones(1,N);   % convert to list if not already
end

Rn=Dn/2;           % need radius and diameter
Ds=min(Dn);        % diameter of smallest particle
Dl=max(Dn);        % diameter of largest particle
if(setDensity)
  L0=sqrt(pi*sum(Dn.^2)/4/phi_set/Gb)/Dl;
else
  L0=L0_set;
end
Ly=L0*Dl;
Lx=Gb*Ly;
  
%% Save Variables
Ek=zeros(Nt,1);     % kinetic energy of particles
Ep=zeros(Nt,1);     % particle-particle potential energy
Ps=zeros(Nt,2,2);    % pressure tensor

if(saveall)
  xs=zeros(Nt,N);     % all particle x-positions
  ys=zeros(Nt,N);     % all particle y-positions
  vxs=zeros(Nt,N);    % all particle x-velocities
  vys=zeros(Nt,N);    % all particle y-velocities
end

%% Initial Conditions
% set up initial velocities with kinetic energy KE0
vx=randn(1,N);               % normal (Maxwell) distribution
vy=randn(1,N);             

vx=vx-mean(vx);              % remove horizontal center of mass motion
vy=vy-mean(vy);              % remove vertical center of mass motion

v2=sum(vx.^2+vy.^2);         % sum of velocities squared
vx=vx.*sqrt(2*KE0./(M*v2));  % rescale so sum(M.*(vx.^2+vy.^2))=2*KE0
vy=vy.*sqrt(2*KE0./(M*v2));

ax_old=0*vx;      % initial particle accelerations for velocity verlet
ay_old=0*vy;

vytp=0;           % initial top plate velocity
aytp_old=0;       % initial top plate accelerations for velocity verlet

% set up particle positions
switch ic
  case 'square'
    Dl=sqrt(Lx*Ly/sum(Dn.^2))*Dl;
    [x y]=ndgrid(Dl/2:Dl:Lx-Dl/2,Dl/2:Dl:Ly-Dl/2); % create grid from N,Dl,Lx,Ly
    ii=randperm(numel(x),N);                       % fill in randomly
    x=reshape(x(ii),1,N);                            
    y=reshape(y(ii),1,N);                            
  case 'hexagonal'
    Dl=min(1.1*Dl,sqrt(.9*4*Lx*Ly/sum(Dn.^2)/pi)*Dl);
    [x y]=ndgrid(Dl/2:Dl:Lx-Dl/2,Dl/2:Dl:2*Ly); % create grid from N,Dl,Lx
    xrange=max(x(:))-min(x(:));         % total range of x
    x=x/(xrange+Dl)*Lx;                 % scale to fit exactly
    [c r]=size(x);                      % number of rows and columns
    if(c>1)
      sx=x(2,1)-x(1,1);                 % x-spacing
    else
      sx=Lx;                            % if only one column the spacing is Lx
    end
    x(:,1:2:end)=x(:,1:2:end)+sx/2;     % shift everyother row
    rs=sqrt(Dl^2-sx.^2/4);              % row shift depends on x-spacing
    y=(y-Dl/2)*rs/Dl+Dl/2;              % shift down by rs per row
    ii=y<=Ly-(sqrt(3)-1)/2*Dl;          % remove ones outside of box
    x=x(ii);
    y=y(ii);
    ii=randperm(numel(x),N);            % fill in randomly
    x=reshape(x(ii),1,N);                            
    y=reshape(y(ii),1,N);                            
  case 'load'
    load(icfile.name,icfile.var);       % load vars from file
end
Dl=max(Dn);                             % set Dl back

%% Setup Plotting
if(plotit)
  clf;                          % clear figure
  xp=mod(x,Lx);                 % put particles in main box 
  yp=mod(y,Ly);                 % nb:periodic images not shown
  h=plotNCirc(xp,yp,Dn,'none'); % plot particles no color
  axis('equal');                % square pixels
  axis([0 Lx 0 Ly]);            % ploting range
  if(pauseit); pause; end       % pause if pauseit true
end

%% Main Loop
Nplotskip=Nplotskip+1;  % fix to skip zero
for nt=1:Nt
  
  % plot particles
  if(rem(nt-1,Nplotskip)==0 && plotit) % plot if nt-1 divides Nplotskip
    xp=mod(x,Lx);                 % put particles in main box
    yp=mod(y,Ly);                 % nb:periodic images not shown
    for np=1:N                         % loop all particles
      % paricles in main box
      set(h(np),'Position',[xp(np)-Rn(np) yp(np)-Rn(np) Dn(np) Dn(np)]);
    end
    drawnow;                           % draw
  end
  
  if(mod(nt,fix(Nt/10))==0)            % visual feedback dot every 10%
    fprintf('.');
  end
  
  % first step in Velocity Verlet integration
  x=x+vx*dt+ax_old.*dt.^2/2;                     % update x
  y=y+vy*dt+ay_old.*dt.^2/2;                     % update y
  
  % save position based quantities after first step
  if (saveall)           % save particle posiitons
    xs(nt,:)=x;
    ys(nt,:)=y;
  end
  
  % Interaction detector and Force Law
  Fx=zeros(1,N);                      % zero forces
  Fy=zeros(1,N);
  
  for nn=1:N                          % check particle n 
    for mm=nn+1:N                     % againt particle m from n+1->N
      dy=y(mm)-y(nn);                 % y-comp of dnm vector from from n->m
      dy=dy-round(dy/Ly)*Ly;          % closest periodic image
      Dnm=Rn(mm)+Rn(nn);              % contact distance
      if(abs(dy)<Dnm)                 % y-distance close enough?    
        dx=x(mm)-x(nn);               % x-comp of dnm vector from from n->m   
        dx=dx-round(dx/Lx)*Lx;        % closest periodic image
        dnm=dx.^2+dy.^2;              % squared distance between n and m
        if(dnm<Dnm^2)                 % overlaping?
          dnm=sqrt(dnm);              % distance
          F=-K*(Dnm/dnm-1);           % force magnitude 
          Fdx=F.*dx;                  % force components
          Fdy=F.*dy;
          Fx(nn)=Fx(nn)+Fdx;           % accumulate x force on n
          Fx(mm)=Fx(mm)-Fdx;           % x force on m (equal opposite) 
          Fy(nn)=Fy(nn)+Fdy;           % accumulate y force on n
          Fy(mm)=Fy(mm)-Fdy;           % y force on m (equal opposite)
          Ep(nt)=Ep(nt)+(Dnm-dnm).^2;  % particle particle potential calc I
          Ps(nt,1,1)=Ps(nt,1,1)+Fdx*dx;  % pressure xx tensot I
          Ps(nt,1,2)=Ps(nt,1,2)+Fdx*dy;  % pressure xy tensot I
          Ps(nt,2,2)=Ps(nt,2,2)+Fdy*dy;  % pressure yy tensot I
        end
      end
    end
  end
  Ep(nt)=K/2*Ep(nt);                  % particle particle potential calc II 
    
  % acceleration with correction for velocity dependent force B
  ax=(Fx./M-B*(vx+ax_old*dt/2))/(1+B*dt/2);     % x-acceleration w/damping
  ay=(Fy./M-B*(vy+ay_old*dt/2))/(1+B*dt/2);   % y-accel. w/damping and g
  
  % second step in Verlet integration
  vx=vx+(ax_old+ax).*dt/2;            % update particle vx  
  vy=vy+(ay_old+ay).*dt/2;            % update particle vy
  
  % save velocity dependent quantities after second setup
  if (saveall)                        % save particle velocity
    vxs(nt,:)=vx;
    vys(nt,:)=vy;
  end
  Ek(nt)=sum(M.*(vx.^2+vy.^2))/2;     % save particle kinetic energy
  
  % save for next step
  ax_old=ax;             
  ay_old=ay; 
end 
%% Final calculations
Ps(:,2,1)=Ps(:,1,2);
Ps=Ps/Lx/Ly;

%% clean up
if (clean)
  fprintf('Done.\n');
end
