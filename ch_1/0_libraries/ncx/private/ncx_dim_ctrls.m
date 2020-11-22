function ncx_dim_ctrls(theHandle,theFrame)

if nargin < 2
  theFrame = 1;
end
fileHandleTag    = ['ncx_files',num2str(theFrame)];
varnameHandleTag = ['ncx_varnames',num2str(theFrame)];
frameHandleTag   = ['ncx_frame',num2str(theFrame)];
dim_ctrlTag      = ['dim_ctrl',num2str(theFrame)];
dim_nameTag      = ['dim_name',num2str(theFrame)];

% get file and varname:
fileHandle = findobj(gcf,'tag',fileHandleTag);
files = get(fileHandle,'userdata');
filev = get(fileHandle,'value');
if iscell(files)
  file = files{filev};
else
  file = files;
end

varnameHandle = findobj(gcf,'tag',varnameHandleTag);
varnames = get(varnameHandle,'string');
varnamev = get(varnameHandle,'value');
if iscell(varnames)
  varname = varnames{varnamev};
else
  varname = varnames;
end

% get reference position of the current frame:
frameHandle = findobj(gcf,'tag',frameHandleTag);
framePos = get(frameHandle,'position');

lengths = getappdata(gcf,'lengths');
w = lengths(1);
h = lengths(2);
d = lengths(3);

x = framePos(1)+d;
y = framePos(2)+framePos(4)-d-h-d-h;


% delete previous:
objs = findobj(gcf,'tag',dim_ctrlTag);
delete(objs);
objs = findobj(gcf,'tag',dim_nameTag);
delete(objs);

if isequal(varname,'-- select --') | isequal(varname,'no var')
  return
end

% find dims:
dim=n_vardim(file,varname);
dimNames = dim.name;
dimLens  = dim.length;

% create uicontrols:
for i=1:length(dimNames)
  show_ctrls(x,y,w,h,dimNames{i},dimLens{i},file,varname,dim_nameTag,dim_ctrlTag,theFrame)
  y = y-h;
end

% apply theme:
ncx_theme('text');
ncx_theme('edit');
if theFrame == 1
  ncx_theme('pushbutton');
  ncx_theme('chind');
end

function show_ctrls(x,y,w,h,dimName,dimLen,file,varname,dim_nameTag,dim_ctrlTag,theFrame)
dim  = uicontrol('tag',dim_nameTag,'style','text','units','normalized','position',[x         y 2*w h],'string',dimName,'HorizontalAlignment','left');
ind1 = uicontrol('tag',dim_ctrlTag,'style','edit','units','normalized','position',[x+2*w     y w   h],'string','1');
ind2 = uicontrol('tag',dim_ctrlTag,'style','edit','units','normalized','position',[x+3*w     y w   h],'string','1');
ind3 = uicontrol('tag',dim_ctrlTag,'style','edit','units','normalized','position',[x+4*w     y w   h],'string',num2str(dimLen));
if theFrame == 1
  dec  = uicontrol('tag',dim_ctrlTag,             'units','normalized','position',[x+5*w     y h   h],'string','<');
  inc  = uicontrol('tag',dim_ctrlTag,             'units','normalized','position',[x+5*w+h   y h   h],'string','>');
  imax = uicontrol('tag',dim_ctrlTag,             'units','normalized','position',[x+5*w+2*h y w   h],'string',num2str(dimLen));

  set(dec, 'userdata',[dim,ind1,ind2,ind3,imax],'callback','ncx(''ncx_chindex(gcbo,''''dec''''   )'')');
  set(inc, 'userdata',[dim,ind1,ind2,ind3,imax],'callback','ncx(''ncx_chindex(gcbo,''''inc''''   )'')');
  set(ind1,'userdata',[dim,ind1,ind2,ind3,imax],'callback','ncx(''ncx_chindex(gcbo,''''manual'''')'')');
  set(ind2,'userdata',[dim,ind1,ind2,ind3,imax],'callback','ncx(''ncx_chindex(gcbo,''''manual'''')'')');
  set(ind3,'userdata',[dim,ind1,ind2,ind3,imax],'callback','ncx(''ncx_chindex(gcbo,''''manual'''')'')');

  set(imax,'callback',['ncx(''ncx_update_dimlen(gcbo,''''',file,''''',''''',varname,''''',''''',dimName,''''')'')']);

  set(dim,'userdata',[ind1,ind2,ind3,imax]);
else
  set(dim,'userdata',[ind1,ind2,ind3]);
end
