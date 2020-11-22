function thedim = n_vardim(file,varname,dimname)
%N_VARDIM   Get dimensions of NetCDF variable
%   Returns the dimension(s) (name and length) of a NetCDF variable.
%
%   Syntax:
%      DIM = N_VARDIM(FILE,VARNAME,DIMNAME)
%
%   Inputs:
%      FILE      NetCDF file
%      VARNAME   Variable
%      DIMNAME   Dimension name [ none ]
%
%   Output:
%      DIM   NetCDF variable dimension length if DINNAME is
%            specified.
%            If not, DIM will be a structure with all the variable
%            dimensions name and associated length.
%            DIM will be empty if the dimension DIMNAME is not found,
%            if the variable VARNAME is not found or if there is an
%            error in the file.
%            Without output and DIMNAME arguments the result is
%            printed.
%
%   Requirement:
%      NetCDF interface for Matlab
%
%   Examples:
%      dims = n_vardim('file.nc','var');
%      dim  = n_vardim('file.nc','temp','time');
%
%   MMA 7-2-2005, martinho@fis.ua.pt
%
%   See also N_VARDIMEXIST, N_FILEDIM

%   Department of Physics
%   University of Aveiro, Portugal

thedim=[];

if nargin < 2
  disp('# arguments required');
  return
end

ncquiet;
nc=netcdf(file);

if isempty(nc)
  return
end

% check varname
if n_varexist(file,varname)
  d=dim(nc{varname});
else
  disp(['# variable ',varname,' not found']);
  close(nc);
  return
end

if isempty(d)
  thedim.name=[];
  thedim.length=[];
else
  for i=1:length(d)
    thedim.name{i}   = name(d{i}); % NetCDF_Dimension
    thedim.length{i} = d{i}(:);    % itsLength
  end
end

close(nc);

if nargin == 3 % look for dimension DIMNAME:
  found = [];
  for i=1:length(thedim.name)
    if isequal(thedim.name{i},dimname)
      found = thedim.length{i};
      break
    end
  end

  thedim = found;
  if isempty(thedim)
    disp(['# dimension ',dimname,' not found in variable ',varname]);
  end
end

% show values if no output argument:
if nargout == 0 & nargin < 3
  fprintf('\n Dimensions in variable %s from NetCDF file\n %s\n\n',varname,file);
  strn = '';
  strl = '';
  for i=1:length(thedim.name)
    strn = strvcat(strn,thedim.name{i});
    strl = strvcat(strl,num2str(thedim.length{i}));
  end
  maxn = size(strn,2);
  maxl = size(strl,2);
  for i=1:length(thedim.name)
    format = ['  --> %',num2str(maxn),'s%',num2str(2*maxl+5-length(num2str(thedim.length{i}))),'s\n'];
    fprintf(1,format,strn(i,:),strl(i,:));
  end
  fprintf(1,'\n');
end
