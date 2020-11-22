function [scale,offset]=n_varscale(file,varname)
%N_VARSCALE   Get scale and offset of NetCDF variable
%   Returns scale_factor and add_offset of var: var=scale*var+offset.
%
%   Syntax:
%      [SCALE,OFFSET]=N_VARSCALE(FILE,VARNAME)
%
%   Inputs:
%      FILE      NetCDF file
%      VARNAME   Variable
%
%   Outputs:
%      SCALE    value of the attribute scale_factor
%      OFFSET   value of the attribute add_offset
%
%   requires;
%      NetCDF interface for Matlab
%
%   Comment:
%      Will try to get scale_factor and add_offset.
%      By default scale = 1 and offset = 0.
%      The output will be empty if there is an error in the file or
%      if the variable does not exists.
%
%   MMA 29-12-2003, martinho@fis.ua.pt
%
%   See also N_VARATT, USE, SHOW

%   Department of Physics
%   University of Aveiro, Portugal

%   07-02-2005 - Improved

scale  = [];
offset = [];

if nargin == 0
   disp('arguments required');
  return
end

% check file:
ncquiet;
nc=netcdf(file);
if isempty(nc)
  return
end
close(nc);

% check varname
if ~n_varexist(file,varname)
  disp(['# variable ',varname,' not found']);
  return
end

if n_varattexist(file,varname,'scale_factor')
  scale  = n_varatt(file,varname,'scale_factor');
else
  scale = 1;
end

if n_varattexist(file,varname,'add_offset')
  offset  = n_varatt(file,varname,'add_offset');
else
  offset = 0;
end
