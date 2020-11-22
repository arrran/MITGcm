function [v,errorstr] = use(file,varname,varargin)
%USE   Load NetCDF variable
%
%   Syntax:
%      [V,ERSTR] = USE(FILE,VARNAME,VARARGIN)
%
%   Inputs:
%      FILE      NetCDF file [ <none> ]
%      VARNAME   NetCDF variable [ <none> ]
%      VARARGIN:
%         DimName, DimValue, Pairs of dimension names and values. The
%            names can be partially when the + sign is at the end or
%            beginning. Ex: for the dimension ocean_time DimName can
%            be '+time'.
%
%   Outputs:
%      V       The data, scale and offset are applied
%      ERSTR   Error string (is displayed if nargin is 1)
%
%   Comments:
%      Works with the NetCDF toolbox available at
%      http://crusty.er.usgs.gov/~cdenham/MexCDF/nc4ml5.html
%      When using the + sign in DimName you may be repeating some
%      indice. Ex: suppose you have the variable dimensions
%      x_rho and eta_rho, if you use '+rho' both the dimensions
%      will have the value DimValue.
%
%   Examples:
%     h = use('grid.nc','h');
%     v = use('v','+time',10,'eta+',50);
%
%   MMA 2001, martinho@fis.ua.pt
%
%   See also SHOW

%   Department of Physics
%   University of Aveiro, Portugal

%   07-02-2005 - Improved
%   23-09-2005 - Added VARARGIN

v        = [];
errorstr = [];

if nargin == 0
  ext = {'*.nc';'*.cdl'};
  [filename, pathname]=uigetfile(ext, 'Choose the NetCDF file');
  file=[pathname,filename];
  if isequal(filename,0)|isequal(pathname,0)
    return
  end
end

if nargin < 2
  varname = input(':: variable ? ','s');
end

vin = varargin;
l=length(vin);
if mod(l,2)
  disp(':: Invalid VARARGIN, use PropertyName,PropertyValue');
  return
end

propNames  = vin(1:2:end-1);
propValues = vin(2:2:end);

if ~n_varexist(file,varname)
  errorstr=sprintf(':: variable %s not found\n',varname);
  disp(errorstr)
  return
end

if nargin > 2
  dims = n_vardim(file,varname);
  if isempty(dims)
    errorstr=sprintf(':: variable %s has unknown dimension',varname);
    return
  end
  dimNames = dims.name;
  dimLen   = dims.length;
  if isempty(dimNames) | isempty(dimLen)
    errorstr=sprintf(':: variable %s has no dimensions',varname);
    return
  end

  for i=1:length(dimNames)
    theValue = ':';
    for j=1:length(propNames)
      if propNames{j}(end) == '+'
        try
          if dimNames{i}(1:length(propNames{j})-1) == propNames{j}(1:end-1)
            theValue = propValues{j};
          end
        end
      elseif propNames{j}(1) == '+'
        try
          if dimNames{i}(end-length(propNames{j})+2:end) == propNames{j}(2:end)
            theValue = propValues{j};
          end
        end
      else
        if isequal(dimNames{i},propNames{j})
          theValue = propValues{j};
        end
      end
    end
    c{i} = theValue;
  end

  str = [];
  for i=1:length(c)
    str = [str,',',num2str(c{i})];
  end
  str=str(2:end);
else
  str = ':';
end

nc = netcdf(file);
  eval(['v = nc{varname}(',str,');'],'v=[];');
nc=close(nc);
v=squeeze(v);

[scale,offset] = n_varscale(file,varname);
v = v*scale + offset;
%if n_varattexist(file,varname,'scale_factor'),  disp(['# scale  of ',varname,' = ',num2str(scale)]);  end
%if n_varattexist(file,varname,'add_offset'),    disp(['# offset of ',varname,' = ',num2str(offset)]); end
