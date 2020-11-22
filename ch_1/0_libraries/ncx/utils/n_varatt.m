function attrib = n_varatt(file,varname,attname)
%N_VARATT   Get attributes of NetCDF variable
%   Recommended standard attributes for variables:
%   (from the NetCDF User's Guide)
%   - units
%   - long_name
%   - valid_min
%   - valid_max
%   - valid range
%   - scale factor
%   - add_offset
%   - _FillValue
%   - missing_value
%   - signedness
%   - FORTRAN_format
%
%   Syntax:
%      ATTRIB = N_VARATT(FILE,VARNAME,ATTNAME)
%
%   Inputs:
%      FILE      NetCDF file
%      VARNAME   variable
%      ATTNAME   attribute [<none>]
%
%   Output:
%     ATTRIB   NetCDF variable attribute value if ATTNAME is
%              specified.
%              If not, ATTRIB will be a structure with all the
%              variable attributes name and value.
%              ATTRIB will be empty if the attribute ATTNAME is not
%              found, if the variable VARNAME is not found or if there
%              is an error in the file.
%              Without output and ATTNAME arguments the result is
%              printed.
%
%   Requirement:
%      NetCDF interface for Matlab
%
%   Examples:
%      attribs = n_varatt('file.nc','varname');
%      attrib  = n_varatt('file.nc','varname','units')
%
%   MMA 25-6-2004, martinho@fis.ua.pt
%
%   See also N_VARATTEXIST, N_FILEATT

%   Department of Physics
%   University of Aveiro, Portugal

%   07-02-2005 - Improved

attrib=[];

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
  a=att(nc{varname});
else
  disp(['# variable ',varname,' not found']);
  close(nc);
  return
end

for i=1:length(a)
    b=a{i};
    attrib.name{i}  = name(b);
    attrib.value{i} = b(:);
end

close(nc);

if nargin == 3 % look for attributen ATTNAME:
  found = [];
  for i=1:length(attrib.name)
    if isequal(attrib.name{i},attname)
      found = attrib.value{i};
      break
    end
  end

  attrib = found;
  if isempty(attrib)
    disp(['# attribute ',attname,' not found in variable ',varname]);
  end
end

% show values if no output argument:
if nargout == 0 & nargin < 3
  fprintf('\n Attributes in variable %s from NetCDF file\n %s\n\n',varname,file);
  strn = '';
  strl = '';
  for i=1:length(attrib.name)
    strn = strvcat(strn,attrib.name{i});

    % deal with numeric values:
    if isnumeric(attrib.value{i})
      if any(size(attrib.value{i}) == 1) & length(attrib.value{i}) == 1
        attrib.value{i} = ['<numeric: ',num2str(attrib.value{i}),'>'];
      else
        attrib.value{i} = ['<numeric: size = ',num2str(size(attrib.value{i})),'>'];
      end

    % deal with long strings:
    else
      if length(attrib.value{i}) > 50;
        attrib.value{i} = [attrib.value{i}(1:50),' ...etc'];
      end
    end

    % deal with empty values:
    if isempty(attrib.value{i})
      attrib.value{i} = '<empty>';
    end

    strl = strvcat(strl,num2str(attrib.value{i}));
  end

  maxn = size(strn,2);
  maxl = size(strl,2);
  for i=1:length(attrib.name)
    format = ['  --> %',num2str(maxn),'s     %s\n'];
    fprintf(1,format,strn(i,:),strl(i,:));
  end
  fprintf(1,'\n');
end
