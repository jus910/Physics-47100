function [h,ht]=plotNCirc(x,y,D,cc,isnum,isclr,prop)
% plotNCirc <Plot N circles>
% Usage:: [h,ht]=plotNCirc(x,y,D,cc[jet(N)],isnum[false],isclr[false],prop[{prop1,val1...}]
%

% revision history:
% 01/06/17 Mark D. Shattuck <mds> plotNCirc.m

%% Parse Input

N=length(x);

iscc=true;
if(~exist('cc','var') || isempty(cc))
	cc=jet(N);
else
  if(ischar(cc))
    if(strcmpi(cc,'none'))
      iscc=false;
    end
  end
end
    

noprop=false;
if(~exist('prop','var') || isempty(prop))
	noprop=true;
end

if(~exist('isnum','var') || isempty(isnum))
	isnum=false;
end

if(~exist('isclr','var') || isempty(isclr))
	isclr=false;
end

if(numel(D)==1)
  D=D*ones(1,N);
end


%% Main

if(isclr)
  clf;
end
  

h=zeros(1,N);
for np=1:N
  h(np)=rectangle(...
    'Position',[x(np)-.5*D(np) y(np)-.5*D(np) D(np) D(np)],...
    'Curvature',[1 1],...
    'edgecolor','b');
  if(iscc)
    set(h(np),'facecolor',cc(np,:));
  end
  if(~noprop)
    set(h(np),prop);
  end
end
axis('equal');
Lx1=min(x-D/2);
Lx2=max(x+D/2);
Ly1=min(y-D/2);
Ly2=max(y+D/2);
axis([Lx1 Lx2 Ly1 Ly2]);

if(isnum)
  ht=text(x,y,num2str((1:N)'));
  set(ht,'horizontalAlign','center');
end




