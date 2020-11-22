function vars = n_filevars(file)
%N_FILEVARS   Get variable names of NetCDF file
%
%   Syntax:
%      VARS = N_FILEVARS(FILE)
%
%   Input:
%      FILE   NetCDF file
%
%   Output:
%      VARS   Variable names
%
%   Example:
%      names = n_filevars('file.nc')
%
%   MMA 24-11-2005, martinho@fis.ua.pt
%
%   See also N_VAREXIST

vars = [];
nc=netcdf(file);
if isempty(nc)
  return
end
v=var(nc);
for i=1:length(v)
  vars{i} = name(v{i});
end
close(nc);

if nargout == 0
  fprintf(1,':: variables in NetCDF file %s\n',file)
  for i=1:length(vars)
    fprintf(1,'  %3g    %s\n',i,vars{i})
 end
end
