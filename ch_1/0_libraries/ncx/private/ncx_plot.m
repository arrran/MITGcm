function ncx_plot

set(gcf,'Pointer','watch');

u=[];
v=[];
w=[];
x=[];
y=[];

% get file of first main var:
theFrame = 1;
[file,varname,theRange,ndims] = get_file_var_range(theFrame); fileU = file; varnameU = varname;
if ndims > 2
  u=[];
elseif ndims~=-1
  % load first main var:
  eval(['u = use(file,varname',theRange,');'],'u=[];');
end
ndimsMain = ndims;
if isempty(u)
  set(gcf,'Pointer','arrow');
  return
end

% get file of 2nd main var:
theFrame = 2;
[file,varname,theRange,ndims] = get_file_var_range(theFrame);
if ndims > 2
  v=[];
elseif ndims~=-1
  % load 2nd main var:
  eval(['v = use(file,varname',theRange,');'],'v=[];');
end
ndimsSec = ndims;

% get file of 3rd main var:
theFrame = 3;
[file,varname,theRange,ndims] = get_file_var_range(theFrame);
if ndims > 2
  w=[];
elseif ndims~=-1
  % load 3rd main var:
  eval(['w = use(file,varname',theRange,');'],'w=[];');
end
ndimsTer = ndims;

% get file of var X:
theFrame = 4;
[file,varname,theRange,ndims] = get_file_var_range(theFrame);
if ndims > 2
  x=[];
elseif ndims~=-1
  % load var X:
  eval(['x = use(file,varname',theRange,');'],'x=[];');
end
ndimsX = ndims;

% get file of var Y:
theFrame = 5;
[file,varname,theRange,ndims] = get_file_var_range(theFrame);
if ndims > 2
  y=[];
elseif ndims~=-1
  % load var Y:
  eval(['y = use(file,varname',theRange,');'],'y=[];');
end
ndimsY=ndims;

% gen XY if needed:
usingXY=1; % used by add coastline, bellow
if isempty(x) | isempty(y)
  x=1:size(u,2);
  y=1:size(u,1);
  [x,y]=meshgrid(x,y);
  usingXY=0;
elseif ~isequal(size(x),size(u)) | ~isequal(size(y),size(u))
  % meshgrid if needed:
  [x,y]=meshgrid(x,y);
end

% get var border:
xb = var_border(x);
yb = var_border(y);

% --------------------------------------------------------------------
% find zoomC
% --------------------------------------------------------------------
cHandle = findobj(gcf,'tag','ncx_zoomC12');
try
  cHandle=cHandle(1);
  xcz = get(cHandle,'xdata');
  ycz = get(cHandle,'ydata');
catch
  xcz = sum(xlim)/2;
  ycz = sum(ylim)/2;
end
cla

% check if is to draw vfield:
isuv=findobj(gcf,'tag','ncx_2ndisuv'); isuv=get(isuv,'value');
isuv=isuv&~isempty(v);
if isuv
  U = sqrt(u.^2+v.^2);

  % get scale: arrow data units per inch
  theScaleHandle = findobj(gcf,'tag','ncx_uvscale'); theScale = str2num(get(theScaleHandle,'string'));
  if ~isnumber(theScale,1)
    % find a nice scale:
    Umean = mean(mean(U));
    axUnits = get(gca,'units'); set(gca,'units','inches'); axPos = get(gca,'position'); set(gca,'units',axUnits);
    ref = size(u).*[axPos(4) axPos(3)];
    theScale = round(max(ref)/Umean);
    set(theScaleHandle,'string',theScale);
  end
  % convert scale inch to scale in data units:
  % axis([min(min(x)) max(max(x)) min(min(y)) max(max(y))]) % use insstead the zoomed one!
  theScale = inch2data(theScale);


  % check if use color:
  useColor = findobj(gcf,'tag','ncx_uvusecolor'); useColor = get(useColor,'value');
  uvColor  = findobj(gcf,'tag','ncx_uvcolor');    uvColor  = get(uvColor,'foregroundColor');
  if useColor
    vfield(x,y,u*theScale,v*theScale,'color',uvColor);
  else
    vfield(x,y,u*theScale,v*theScale,U);
  end

else
  % ---------------------------- main var:
  % get if pcolor and contour values:
  ispcolor = findobj(gcf,'tag','ncx_1stispcolor'); ispcolor = get(ispcolor,'value');
  if ~ispcolor
    cvals  = findobj(gcf,'tag','ncx_1stcontVals');  cvals  = str2num(get(cvals, 'string'));
    if ~isnumber(cvals)
      cvals=10;
    end
  end

  if ndimsMain == 1
    try
      plot(x,u);
    catch
      plot(u);
    end
  elseif ndimsMain == 2
    if ispcolor
      try
        pc = pcolor(x,y,u);
      catch
        pc = pcolor(u);
      end
      set(pc,'faceColor','interp','edgeColor','none');
    else
      try
        [cs,ch] = contour(x,y,u,cvals);
      catch
        [cs,ch] = contour(u,cvals);
      end
      % store contour handles:
      cont_var1.cs=cs;
      cont_var1.ch=ch;
      setappdata(gcf,'contours_1stvar',cont_var1);
    end
  end

  % --------------------------- 2nd var
  ish = ishold; hold on
  % get contour values and color:
  cvals  = findobj(gcf,'tag','ncx_2ndcvals');  cvals  = str2num(get(cvals, 'string'));
  ccolor = findobj(gcf,'tag','ncx_2ndccolor'); ccolor = get(ccolor,'foregroundColor');
  if ~isnumber(cvals)
    cvals=10;
  end

  if ndimsSec == 1
    try
      plot(x,v);
    catch
      plot(v);
    end
  elseif ndimsSec == 2
    try
      [cs,ch]=contour(x,y,v,cvals,'k'); set(ch,'color',ccolor);
    catch
      [cs,ch]=contour(v,cvals,'k'); set(ch,'color',ccolor);
    end
    % store contour handles:
    cont_var2.cs=cs;
    cont_var2.ch=ch;
    setappdata(gcf,'contours_2ndvar',cont_var2);
  end
  if ~ish, hold off, end

end % isuv

% ----------------------------- 3rd var
ish = ishold; hold on
% get contour values and color:
cvals  = findobj(gcf,'tag','ncx_3rdcvals');  cvals  = str2num(get(cvals, 'string'));
ccolor = findobj(gcf,'tag','ncx_3rdccolor'); ccolor = get(ccolor,'foregroundColor');

% get if use pcolor:
ispcolor = findobj(gcf,'tag','ncx_3rdispcolor'); ispcolor = get(ispcolor,'value');

if ~isnumber(cvals)
  cvals=10;
end

if ndimsTer == 1
  try
    plot(x,w);
  catch
    plot(w)
  end
elseif ndimsTer == 2
  if ispcolor
    try
      pc=pcolor(x,y,w);
    catch
      pc=pcolor(w);
    end
    set(pc,'faceColor','interp','edgeColor','none');
  else
    try
      [cs,ch]=contour(x,y,w,cvals,'k'); set(ch,'color',ccolor);
    catch
      [cs,ch]=contour(w,cvals,'k'); set(ch,'color',ccolor);
    end
    % store contour handles:
    cont_var3.cs=cs;
    cont_var3.ch=ch;
    setappdata(gcf,'contours_3rdvar',cont_var3);
  end
end
if ~ish, hold off, end

% --------------------------------------------------------------------
% plot var border:
% --------------------------------------------------------------------
if ndimsMain > 1
  ish = ishold; hold on
  DataBorder = plot(xb,yb,'b:','tag','ncx_DataBorder');
  if ~ish, hold off; end
end

% --------------------------------------------------------------------
% caxis:
% --------------------------------------------------------------------
theCaxis = findobj(gcf,'tag','ncx_Caxis');  theCaxis = str2num(get(theCaxis,'string'));
try
  caxis(theCaxis);
catch
  caxis('auto');
end

% --------------------------------------------------------------------
% colorbar:
% --------------------------------------------------------------------
% colorbar is used with patches or surfaces:
existPatches  = ~isempty(findobj(gca,'type','patch'));
existSurfaces = ~isempty(findobj(gca,'type','surface'));
cb = findobj(gcf,'tag','ncx_colorbar');
cbc = get(cb,'children');
if existPatches | existSurfaces
  set(cb,'visible','on');
  set(cb,'userdata',[]) % needed for R>12
  colorbar(cb);
  set(cb,'tag','ncx_colorbar');
else
  set(cb,'visible','off');
  delete(cbc);
end

% --------------------------------------------------------------------
% coastline:
% --------------------------------------------------------------------
if usingXY & ndimsMain > 1
  ncx('ncx_coastl(''add'')');
end

% --------------------------------------------------------------------
% deal with changes on PlotBoxAspectRatio when axis equal,
% in the case of isuv:
% --------------------------------------------------------------------
if isuv
  xdb = get(DataBorder,'xdata');
  ydb = get(DataBorder,'ydata');
  axUnits = get(gca,'units');
  set(gca,'units','pixels');
  axpos = get(gca,'position');
  set(gca,'units',axUnits);
  a = max(xdb)-min(xdb);
  b = max(ydb)-min(ydb);
  H = axpos(4);
  W = axpos(3);
  r = (b/a)*(W/H);
  x1 = min(xdb);
  x2 = max(xdb);
  y1 = min(ydb);
  y2 = max(ydb);
  xc = (x1+x2)/2;
  yc = (y1+y2)/2;
  dx = x2-x1;
  dy = y2-y1;
  if r > 1
    k = r;
    x2 = xc+k*dx/2;
    x1 = xc-k*dx/2;
  elseif r < 1
    k = 1/r;
    y2 = yc+k*dy/2;
    y1 = yc-k*dy/2;
  end
  xdb = [xdb nan x1 nan x2];
  ydb = [ydb nan y1 nan y2];
  set(DataBorder,'xdata',xdb,'ydata',ydb);
end

% --------------------------------------------------------------------
% apply zoom, and show zoomC:
% --------------------------------------------------------------------
zoomC = findobj(gcf,'tag','ncx_zoomCentre');
if get(zoomC,'value')
  zoomc('add',[xcz,ycz]);
end
if ndimsMain > 1
  ncx('ncx_zoom');
end
if isuv
  show_arrowScale(U,fileU,varnameU,theScale)
end

% --------------------------------------------------------------------
% apply axes theme:
% --------------------------------------------------------------------
ncx_theme('axes');

set(gcf,'Pointer','arrow');

function [file,varname,theRange,ndims] = get_file_var_range(theFrame)
frameTag = ['ncx_frame',num2str(theFrame)];
theFrameHandle = findobj(gcf,'tag',frameTag);
isVisible = get(theFrameHandle,'visible');

fileHandleTag    = ['ncx_files',num2str(theFrame)];
varnameHandleTag = ['ncx_varnames',num2str(theFrame)];
dim_nameTag      = ['dim_name',num2str(theFrame)];

% get file:
fileHandle = findobj(gcf,'tag',fileHandleTag);
files = get(fileHandle,'userdata');
filev = get(fileHandle,'value');
if iscell(files)
  file = files{filev};
else
  file = files;
end

% get varname:
varnameHandle = findobj(gcf,'tag',varnameHandleTag);
varnames = get(varnameHandle,'string');
varnamev = get(varnameHandle,'value');
if iscell(varnames)
  varname = varnames{varnamev};
else
  varname = varnames;
end

if isequal(isVisible,'off')
  theRange = '';
  ndims = -1;
  return
end

if isequal(varname,'-- select --')
  theRange = '';
  ndims = -1;
  return
end

% get range:
theRange = '';
objs = findobj(gcf,'tag',dim_nameTag);
ndims = 0;
for i=1:length(objs)
  handles = get(objs(i),'userdata');
  dimName = get(objs(i),'string');
  ind1 = handles(1); ind1=get(ind1,'string');
  ind2 = handles(2); ind2=get(ind2,'string');
  ind3 = handles(3); ind3=get(ind3,'string');

  theRange = [theRange,',', sprintf('''%s'',''%s:%s:%s''',dimName,ind1,ind2,ind3)];

  if length([str2num(ind1):str2num(ind2):str2num(ind3)]) > 1
    ndims = ndims+1;
  end
end

function show_arrowScale(U,fileU,varnameU,theScale)
xl = xlim;
yl = ylim;
scx = xl(1)+diff(xl)/3;
scy = yl(1)-diff(yl)/4.5;

hold on
Usc = mean(mean(U));
[Usc,iv,e] = vscale(Usc);
Usc=2*Usc*10^e;
scHandle = vfield(scx,scy,Usc*theScale,0); set(scHandle,'tag','ncx_vfieldScale','Clipping','off');

% find u units
if n_varattexist(fileU,varnameU,'units')
  units = n_varatt(fileU,varnameU,'units');
else
  units = 'units ??';
end
str = ['  ',num2str(Usc),'  (',units,')'];
scText   = text(scx+Usc*theScale,scy,str);
hold off
