% function circle(x,y,r,spec)
function circle(varargin)

% plots a circle.
% x and y are the coordinates of the center of the circle
% r is the radius of the circle
% 0.01 is the angle step, bigger values will draw the circle faster but
% you might notice imperfections (not very smooth)

% cf changed to varargin for z dimension

if nargin==3
    x = varargin{1};
    y = varargin{2};
    z = NaN;
    r = varargin{3};
    spec = 'k-';
elseif nargin==4
    x = varargin{1};
    y = varargin{2};
    z = NaN;
    r = varargin{3};
    spec = varargin{4};
elseif nargin==5
    x = varargin{1};
    y = varargin{2};
    z = varargin{3};
    r = varargin{4};
    spec = varargin{5};
else
    error('incorrect nargin');
end
    
ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);

if nargin==5
    plot3(x+xp,y+yp,ones(size(xp))*z,spec);
else
    plot(x+xp,y+yp,spec);
end


end