function is=isnumber(v,n,space)
%ISNUMBER   Enhanced version of isnumeric
%   Similar to isnumeric but checks length/size, NaNs and range
%   (space).
%
%   Syntax:
%      IS = ISNUMBER(V,N,SPACE)
%      IS = ISNUMBER(V,SPACE)
%
%   Inputs:
%      V       Data, N-D array
%      N       Desired length or size
%      SPACE   Allowed values (any by default):
%              'Z', 'Z+', 'Z0+', 'Z-', 'Z0-', 'R+', 'R0+', 'R-', 'R0-'
%
%   Output:
%      IS   0 or 1
%
%   Examples:
%      isnumber([1 2 3],3)       % 1
%      isnumber([1 2 3])         % 1
%      isnumber('www')           % 0
%      isnumber([1 2.1],'Z')     % 0
%      isnumber(rand(3,4),[3 4]) % 1
%
%   MMA, 8-2003, martinho@fis.ua.pt

%   Department of Physics
%   University of Aveiro, Portugal

%   02-02-2005 - Introduced argument SPACE
%   11-10-2005 - Allow N-D arrays, instead of vectors

is=0;
check_size  = 0;
check_space = 0;

if nargin == 2
  if isstr(n)
    space = n;
    check_space = 1;
  else
    check_size = 1;
  end
end
if nargin ==3
  check_size  = 1;
  check_space = 1;
end

if isempty(v),    return; end
if ~isnumeric(v), return, end
if any(isnan(v)), return, end

if check_size
  if length(n) == 1
    if ~(ndims(v) == 2 & any(size(v)==1) & length(v) == n), return, end
  else
    if ~(isequal(size(v),n)), return, end
  end
end

if check_space
  switch space
    case 'Z',   if ~all( v == floor(v)          ), return, end
    case 'Z-',  if ~all( v == floor(v) & v <  0 ), return, end
    case 'Z0-', if ~all( v == floor(v) & v <= 0 ), return, end
    case 'Z+',  if ~all( v == floor(v) & v >  0 ), return, end
    case 'Z0+', if ~all( v == floor(v) & v >= 0 ), return, end

    case 'R-',  if ~all(v <  0), return, end
    case 'R0-', if ~all(v <= 0), return, end
    case 'R+',  if ~all(v >  0), return, end
    case 'R0+', if ~all(v >= 0), return, end
    otherwise, disp([':: invalid space : ',space]); return
  end
end

is=1;
