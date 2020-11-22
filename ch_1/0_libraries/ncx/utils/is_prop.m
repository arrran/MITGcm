function res=is_prop(handle,prop)
%IS_PROP   Check handle property
%
%   Syntax:
%      RES = IS_PROP(HANDLE,PROP)
%
%   Inputs:
%      HANDLE   Handle to check
%      PROP     Property to find
%
%   Output:
%      RES   0 or 1
%
%   Example:
%      p=plot(1:10);
%      res=is_prop(p,'Color')
%
%   MMA 13-8-2003, martinho@fis.ua.pt

%   Department of physics
%   University of Aveiro

%   2007-12-12 - Renamed from isprop due to conflicts with a MATLAB
%                built in function

res=0;
if ~ishandle(handle)
  return
end
s=get(handle);
if isfield(s,prop)
  res=1;
end
