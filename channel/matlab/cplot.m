function mov = cplot( varargin )
% cplot is to play 2d movie for the dynamics of certain contour and return
% the movie object.
% 
% The input usually involves two matrices, the first one for x, the second 
% one for y, the dimensions of which should match. The program plots the
% x-y figure in a coloumn first manner. In other words, the movie plots a
% certain column of y against the corresponding column of x for a certain
% frame
%
% usage: 'n' or 'new': play the movie in a new frame
%        'l' or 'label': indicate the x-label and y-label in the following
%                        two arguments
%        'p' or 'pause': vary the pause time from default 0.2s in the
%                        following argument
%        'a' or 'avi':   save the movie in avi format
%        'title':        movie title
%        'name':         determine the avi file name in the following
%                        argument
%        '[]':           determine coloumn first(1) or row first(2)
%        'color':        point color on the graph
%
%   creator: Sirius
%   date: 4/18/2012
%   copy right reserved

vinput=0;                   % number of input matries, 1 or 2
fhandle=@plot;              % function handle, can be @plot, @scatter, etc.
xlb='';                     % x-label
ylb='';                     % y-label
pause_time=0.2;             % pause time for the movie
newgraph=0;                 % whether or not the movie is played in a new frame
avi=0;                      % whether or not the movie is saved in avi format
title_name='';              % movie title name
avi_name='';                % avi file name
point_color='k';            % point color on the graph

if nargin>=2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
    x=varargin{1};
    y=varargin{2};
    vinput=2;
elseif nargin>=1 && isnumeric(varargin{1})
    y=varargin{1};
    x=(1:size(y,1))';
    vinput=1;
else
    error('invalid input')
end

if nargin>vinput
    i=vinput+1;
    while(i<=nargin)
        j=i;
        if nargin-i>=2 && strcmpi(varargin{i},'l') || strcmpi(varargin{i},'label')
            xlb=varargin{i+1};
            ylb=varargin{i+2};
            j=i+2;
        elseif isa(varargin{i},'function_handle')
            fhandle=varargin{i};
        elseif strcmpi(varargin{i},'n') || strcmpi(varargin{i},'new')
            newgraph=1;
        elseif strcmpi(varargin{i},'p') || strcmpi(varargin{i},'pause')
            pause_time=varargin{i+1};
            j=i+1;
        elseif strcmpi(varargin{i},'a') || strcmpi(varargin{i},'avi')
            avi=1;
        elseif strcmpi(varargin{i},'name')
            avi_name=varargin{i+1};
            j=i+1;
        elseif strcmpi(varargin{i},'title')
            title_name=varargin{i+1};
            j=i+1; 
        elseif strcmpi(varargin{i},'[]')
            if varargin{i+1}==2
                x=x';
                y=y';
            end
            j=i+1;
        elseif strcmpi(varargin{i},'color')
            point_color=varargin{i+1};
            j=i+1;
        end
        i=j+1;
    end
end

if sum(size(x)==size(y))==0
    error('dimensions of the input matrices mismatch')
elseif sum(size(x)==size(y))==1
    if size(x,1)~=size(y,1)
        x=x';
        y=y';
    end
    
    if size(x,2)==1
        x=x*ones(1,size(y,2));
    elseif size(y,2)==1
        y=y*ones(1,size(x,2));
    else
        error('dimensions of the input matrices mismatch')
    end
end

xmin=min(min(x));
xmax=max(max(x));
ymin=min(min(y));
ymax=max(max(y));

if newgraph~=0
    figure;
else
    clf;
end
frame=size(x,2);

mov(1:frame)=struct('cdata',[],'colormap',[]);

xlabel(xlb);
ylabel(ylb);
title(title_name);
axis([xmin xmax ymin ymax])
set(gca,'NextPlot','replacechildren')

for i=1:frame
    fhandle(x(:,i),y(:,i),point_color);   
    
    mov(i)=getframe;
    pause(pause_time)
end

% save to avi file
if avi~=0
    if isempty(xlb) && isempty(ylb) && isempty(avi_name)
        avi_name=[xlb '_' ylb];
    end
    if exist([avi_name 'avi'],'file')
        avi_name=[avi_name '1'];
    end
    movie2avi(mov, [avi_name '.avi'], 'compression', 'None');
    disp([avi_name '.avi saved!'])
end

end

