function out = n_varexist(file,varname)
%N_VAREXIST   Check if NetCDF variable exists
%
%   Syntax:
%      OUT = N_VAREXIST(FILE,VARNAME)
%
%   Inputs:
%      FILE      NetCDF file
%      VARNAME   Variable
%
%   Output:
%      OUT   1 = yes; 0 = no; [] = error in file
%
%   Requirement:
%      NetCDF interface for Matlab
%
%   Example:
%      out = n_varexist('filename.nc','lon');
%
%   MMA 7-6-2004, martinho@fis.ua.pt

%   Department of Physics
%   University of Aveiro, Portugal

out=[];

ncquiet;
nc=netcdf(file);
if isempty(nc)
  return
end

vars=var(nc);

out=0;
for i=1:length(vars)
  if isequal(name(vars{i}),varname)
    out=1;
    close(nc);
    return
  end
end

close(nc);
