function [lonc,latc] = rot_longitude(lon,lat)
%ROT_LONGITUDE   Swap longitude between -180:180 and 0:360
%   Conversion between -180:180 to 0:360 is done if negative values
%   are found in longitude, and the inverse, otherwise.
%
%   Syntax:
%      [LONC,LATC] = ROT_LONGITUDE(LON,LAT)
%
%   Inputs:
%     LON, LAT   Longitude and latitude vectors
%
%   output:
%      LATC, LONC Modified longitude and latitude (nans are added)
%
%   Example:
%      load coastline % must have lon and lat variables
%      figure, plot(lon,lat); axis equal
%      [lonc,latc] = rot_longitude(lon,lat);
%      figure plot(lonc,latc), axis equal
%
%   MMA 23-10-2004, martinho@fis.ua.pt

%   Department of Physics
%   University of Aveiro, Portugal

lonc = [];
latc = [];

if nargin < 2
  disp('## arguments needed...');
  return
end

is_neg = 0;
if any(lon<0)
  is_neg = 1;
end

if is_neg      % conditions to change from -180:180 to 0:360:
  add  = 360;
  add2 = 0;

  cnd  = lon < 0;
  cnd2 = 'lonc < 0;';
else           % conditions to change from 0:360 to -180:180:
  add  = 360;
  add2 = -360;

  cnd  = lon < 180;
  cnd2 = 'lonc < 180;';
end

d = diff(cnd);
d = find(d); d = reshape(d,1,length(d));
d=[0 d length(lon)];

for i=1:length(d)-1
  lo = lon(d(i)+1:d(i+1)); l = length(lo);
  la = lat(d(i)+1:d(i+1));
  lonc(end+1:end+l) = lo;
  lonc(end+1) = nan;

  latc(end+1:end+l) = la;
  latc(end+1) = nan;
end
lonc = lonc(1:end-1); % last nan is not needed
latc = latc(1:end-1);

eval(['cnd = ',cnd2]);
lonc(cnd) = lonc(cnd) + add;
lonc = lonc + add2;
