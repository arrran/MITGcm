function ncx_load

[filename, pathname] = uigetfile({'*.nc;*.cdf;', 'NetCDF Files (*.nc, *.cdf)'},'Pick a NetCDF file');
if isequal(filename,0) | isequal(pathname,0)
  return
else
  file = fullfile(pathname, filename);
end

obj = findobj(gcf,'tag','ncx_files1');
ud  = get(obj,'userdata');
str = get(obj,'string');

ud{end+1} = file;
set(obj,'userdata',ud);

strMax = 20;
strLen = length(file);
new = file(strLen-min(strMax,strLen)+1:end);
if isempty(str) | isequal(str,'[0] no file loaded yet')
  n=1;
  str={};
else
  n = length(str)+1;
end
new = [['[',num2str(n),']'],' -- ',new];
str{end+1} = new;
set(obj,'string',str);

% after update the main files list, update all the others.
% also update the vars list if is the unique file loaded!!

if 1
  try, handleFiles(1)     = findobj(gcf,'tag','ncx_files2'); catch, handleFiles(1)     = nan; end
  try, handleFiles(end+1) = findobj(gcf,'tag','ncx_files3'); catch, handleFiles(end+1) = nan; end
  try, handleFiles(end+1) = findobj(gcf,'tag','ncx_files4'); catch, handleFiles(end+1) = nan; end
  try, handleFiles(end+1) = findobj(gcf,'tag','ncx_files5'); catch, handleFiles(end+1) = nan; end

  for i=1:n
    filesn{i} = sprintf('[%g]',i);
  end

  for i=1:length(handleFiles)
    if ~isnan(handleFiles(i))
      set(handleFiles(i),'string',filesn);
      set(handleFiles(i),'userdata',ud);
    end
  end
end

if n==1
  try, handleVars(1)     = findobj(gcf,'tag','ncx_varnames1'); catch, handleVars(1)     = nan; end
  try, handleVars(end+1) = findobj(gcf,'tag','ncx_varnames2'); catch, handleVars(end+1) = nan; end
  try, handleVars(end+1) = findobj(gcf,'tag','ncx_varnames3'); catch, handleVars(end+1) = nan; end
  try, handleVars(end+1) = findobj(gcf,'tag','ncx_varnames4'); catch, handleVars(end+1) = nan; end
  try, handleVars(end+1) = findobj(gcf,'tag','ncx_varnames5'); catch, handleVars(end+1) = nan; end

  varnames = n_filevars(file);
  varnames(2:end+1) = varnames;
  varnames{1} = '-- select --';

  for i=1:length(handleVars)
    if ~isnan(handleVars(i))
      set(handleVars(i),'value',1);
      set(handleVars(i),'string',varnames);
    end
  end
end


