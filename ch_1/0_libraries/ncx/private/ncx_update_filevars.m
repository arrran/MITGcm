function ncx_update_filevars(theHandle)


tag = get(theHandle,'tag');
theFrame = tag(length('ncx_files')+1:end);

varnameHandleTag = ['ncx_varnames',theFrame];
varnameHandle    = findobj(gcf,'tag',varnameHandleTag);
dim_ctrlTag      = ['dim_ctrl',num2str(theFrame)];
dim_nameTag      = ['dim_name',num2str(theFrame)];

% get file:
fileHandle = theHandle;
files = get(fileHandle,'userdata');
filev = get(fileHandle,'value');
if iscell(files)
  file = files{filev};
else
  file = files;
end

% get varnames:
varnames = n_filevars(file);
varnames(2:end+1) = varnames;
varnames{1} = '-- select --';

set(varnameHandle,'value',1);
set(varnameHandle,'string',varnames);

% delete current controls:
objs = findobj(gcf,'tag',dim_ctrlTag);
delete(objs);
objs = findobj(gcf,'tag',dim_nameTag);
delete(objs);

