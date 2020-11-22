function out =  n_varattexist(file,varname,attname)
%N_VARATTEXIST   Check if NetCDF variable attribute exists
%
%   Syntax:
%      OUT = N_VARATTEXIST(FILE,VARNAME,ATTNAME)
%
%   Inputs:
%      FILE      NetCDF file
%      VARNAME   Variable
%      ATTNAME   Attribute name
%
%   Output:
%      OUT   1 = yes; 0 = no; [] = error in file or if VARNAME not
%            found
%
%   Requirement:
%      NetCDF interface for Matlab
%
%   Example:
%      out = n_varattexist('file.nc','temp','units')
%
%   MMA 7-2-2005, martinho@fis.ua.pt
%
%   See also N_VARATT, N_FILEATT, N_FILEATTEXIST

%   Department of Physics
%   University of Aveiro, Portugal

out=[];

if nargin < 3
  disp('# arguments required');
  return
end

attrib = n_varatt(file,varname);

if isempty(attrib)
  return
end

out=0;
for i=1:length(attrib.name)
  if isequal(attrib.name{i},attname)
    out = 1;
    return
  end
end
