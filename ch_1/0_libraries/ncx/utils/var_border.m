function [x,xc] = var_border(M)
%VAR_BORDER   Get border of 2D array
%
%   Syntax:
%      [X,XC] = VAR_BORDER(M)
%
%   Input:
%      M   2D array
%
%   Output:
%      X    Border
%      XC   Values at the 4 corners
%
%   Example:
%      x = 1:10;
%      y = 1:10;
%      [x,y] = meshgrid(x,y);
%      [xx,xxc] = var_border(x);
%      [yy,yyc] = var_border(y);
%      M = rand(10,10);
%      [m,mc] = var_border(M);
%      figure
%      plot(xx,yy); hold on
%      plot3(xx,yy,m,'r')
%      view([-30 60])
%      plot3(xxc,yyc,mc,'bo')
%      axis([-1 11 -1 11 -1 2])
%
%   MMA 18-8-2004, martinho@fis.ua.pt
%
%   See also PLOT_BORDER3D, ROMS_BORDER

%   Department of Physics
%   University of Aveiro, Portugal

x  = [];
xc = [];

if nargin == 0
  disp('» no variable')
  return
end


xl = M(:,1);
xt = M(end,:);  xt = xt';
xr = M(:,end);  xr = flipud(xr);
xb = M(1,:);    xb = flipud(xb');

x =  [xl; xt; xr; xb];

% corners:
xc =  [xl(1) xl(end) xr(1) xr(end)];
