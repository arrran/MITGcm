function ncx(varargin)
%NCX   NetCDF eXplorer
%   NCX is an interface for NetCDF visualisation.
%   It is an improvement of NCDView with many new features, intuitive
%   usage and with a more clear, easy but faster code.
%   NCX already includes the Matlab/NetCDF interface for Matlab
%   versions from 5.3 until the current release 14, for Windows and
%   Linux (authors: Charles R. Denham and John Evans).
%
%   For a tuturial see the site:
%      http://neptuno.fis.ua.pt/~mma/ncx
%
%   If you find problems or have doubts/comments, use my email:
%      martinho@fis.ua.pt
%
%   MMA 12-2005

%   Department of Physics
%   University of Aveiro, Portugal

if isempty(varargin)

% ------------------------------------------------------------------
% figure...
% ------------------------------------------------------------------
  prevNCX = findobj(0,'name','NCX');
  close(prevNCX);

  screenSize  = get(0,'ScreenSize');
  screenUnits = get(0,'Units');
  set(0,'Units','pixels');
  screenSize  = get(0,'ScreenSize');
  set(0,'Units',screenUnits);
  screenCx = screenSize(3)/2;
  screenCy = screenSize(4)/2;
  guiW = 800;
  guiH = 530;
  guiX = screenCx-guiW/2;
  guiY = screenCy-guiH/2;
  figure('units','pixel','NumberTitle','off','Name','NCX','MenuBar','none','position',[guiX guiY guiW guiH],'CloseRequestFcn','ncx(''ncx_close'')');

  x = 0.01; x0=x;
  y = 0.99;
  w = 0.045;
  h = 0.0315;
  d = 0.01;

  lengths = [w h d];
  setappdata(gcf,'lengths',lengths); % needed in ncx_dim_ctrls.m

% ------------------------------------------------------------------
% set path
% ------------------------------------------------------------------
  ncx_paths

% --------------------------------------------------------------------
% menus:
% --------------------------------------------------------------------
  % menus:
  menuFile     = uimenu('label','file');
  menuFileLoad = uimenu('parent',menuFile,'label','load','callback','ncx(''ncx_load'')');
  menuAddCoast = uimenu('parent',menuFile,'label','add coastline','callback','ncx(''ncx_coastl(gcbo)'')','checked','off','tag','ncx_addCoastline','checked','on','separator','on');

  menuNcx        = uimenu('label','ncx');
  menuNcxBasic   = uimenu('parent',menuNcx,'label','basic',      'callback','ncx(''ncx_modes(gcbo,1)'')','tag','ncx_menu_modes');
  menuNcxBasicXy = uimenu('parent',menuNcx,'label','basic+xy',   'callback','ncx(''ncx_modes(gcbo,2)'')','tag','ncx_menu_modes','checked','on');
  menuNcxTwo     = uimenu('parent',menuNcx,'label','two vars',   'callback','ncx(''ncx_modes(gcbo,3)'')','tag','ncx_menu_modes');
  menuNcxThree   = uimenu('parent',menuNcx,'label','three vars', 'callback','ncx(''ncx_modes(gcbo,4)'')','tag','ncx_menu_modes');

  menuAbout      = uimenu('parent',menuNcx,'label','about', 'callback','ncx(''ncx_about'')','separator','on');

% --------------------------------------------------------------------
% frame 1:
% --------------------------------------------------------------------
  H = 7*h+2*d + d + d;
  W = 6*w+2*h + 2*d;   Wmax = W;

  frameDims  = uicontrol('style','frame',    'units','normalized','position', [x       y-H W H],'tag','ncx_frame1');
  fileHandle = uicontrol('style','popupmenu','units','normalized','position', [x+d     y-d-h       5*w     h      ],'string','[0] no file loaded yet','tag','ncx_files1',...
                         'callback','ncx(''ncx_update_filevars(gcbo)'')');
  varHandle  = uicontrol('style','popupmenu','units','normalized','position', [x+5*w+d y-d-h       2*h+w   h      ],'string','no var','tag','ncx_varnames1',...
    'callback',['ncx(''ncx_dim_ctrls(gcbo)'')']);

  dsp   = uicontrol('units','normalized','position', [x+W-w-d-1*h y-H+d w+1*h h],'string','disp');
  dspcb = uicontrol('style','checkbox','units','normalized','position', [x+W-1.5*w-d-1*h y-H+d w/2 h],'string','','tag','ncx_dispcb');
  set(dsp,'userdata',dspcb,'callback','ncx(''ncx_plot'')');

  is1pcolor   = uicontrol('style','radiobutton','units','normalized','position', [x+d y-H+d+h 2*w h],'string','pcolor','tag','ncx_1stispcolor');
  is1cont     = uicontrol('style','radiobutton','units','normalized','position', [x+d y-H+d 2*w h],'string','contour','tag','ncx_1stiscont');
  is1contVals = uicontrol('style','edit','units','normalized','position', [x+d+2*w y-H+d 2*w h],'string','contours','tag','ncx_1stcontVals');
  set(is1pcolor, 'callback','ncx(''rbuttons(gcbo,{''''ncx_1stispcolor'''',''''ncx_1stiscont''''})'')','value',1);
  set(is1cont,   'callback','ncx(''rbuttons(gcbo,{''''ncx_1stispcolor'''',''''ncx_1stiscont''''})'')','value',0);
  cont1clab   = uicontrol(              'units','normalized','position', [x+d+4*w y-H+d w h],'string','clab','callback','ncx(''ncx_contlabel(gcbo,1)'')');


% --------------------------------------------------------------------
% frame 2:
% --------------------------------------------------------------------
  y = y-H - d;
  H = 6*h+2*d + d;
  W = 5*w+2*d;
  frame2     = uicontrol('style','frame',    'units','normalized','position', [x       y-H      W H],'tag','ncx_frame2');
  theFiles2  = uicontrol('style','popupmenu','units','normalized','position', [x+d     y-d-h       2*w     h      ],'string','[0]','tag','ncx_files2',...
                        'callback','ncx(''ncx_update_filevars(gcbo)'')');
  theVars2   = uicontrol('style','popupmenu','units','normalized','position', [x+2*w+d y-d-h       2*h+w h          ],'string','no var','tag','ncx_varnames2',...
                        'callback',['ncx(''ncx_dim_ctrls(gcbo,2)'')']);

  frame22    = uicontrol('style','frame',    'units','normalized','position', [x+W-d/5       y-H     Wmax-W+d/5 H],'tag','ncx_frame22');
  isuv       = uicontrol('style','radiobutton','units','normalized','position', [x+W+d       y-h-d     2*w h],'string','uv field','tag','ncx_2ndisuv');
  iscont     = uicontrol('style','radiobutton','units','normalized','position', [x+W+d       y-h-d-h   2*w h],'string','contour','tag','ncx_2ndiscont');
  set(isuv,  'callback','ncx(''rbuttons(gcbo,{''''ncx_2ndisuv'''',''''ncx_2ndiscont''''})'')','value',0);
  set(iscont,'callback','ncx(''rbuttons(gcbo,{''''ncx_2ndisuv'''',''''ncx_2ndiscont''''})'')','value',1);
  contVals   = uicontrol('style','edit',       'units','normalized','position', [x+W+d       y-h-d-h-h 2*w h],'string','contours','tag','ncx_2ndcvals');
  contcolor  = uicontrol(                      'units','normalized','position', [x+W+d       y-4*h-d   0.5*w h],'string','c','callback','ncx(''color_button(gcbo)'')','tag','ncx_2ndccolor');
  contclab   = uicontrol(                      'units','normalized','position', [x+W+d+w/2   y-4*h-d   1.5*w h],'string','clab','callback','ncx(''ncx_contlabel(gcbo,2)'')');

  uvusecolor = uicontrol('style','checkbox','units','normalized','position', [x+W+d       y-5*h-2*d   0.5*w h],'string','','callback','','tag','ncx_uvusecolor');
  uvcolor    = uicontrol(                      'units','normalized','position', [x+W+d+w/2   y-5*h-2*d   0.5*w h],'string','c','callback','ncx(''color_button(gcbo)'')','tag','ncx_uvcolor');
  uvscale    = uicontrol('style','edit',       'units','normalized','position', [x+W+d+1*w   y-5*h-2*d   1*w h],'string','uv scale','tag','ncx_uvscale');

  frame2Handles = [frame2;theFiles2;theVars2;frame22;isuv;iscont;contVals;contcolor;contclab;uvusecolor;uvcolor;uvscale];


% --------------------------------------------------------------------
% frame 3:
% --------------------------------------------------------------------
  y=y-H - d;
  frame3    = uicontrol('style','frame',    'units','normalized','position', [x       y-H     W H],'tag','ncx_frame3');
  theFiles3 = uicontrol('style','popupmenu','units','normalized','position', [x+d     y-d-h       2*w     h      ],'string','[0]','tag','ncx_files3',...
                        'callback','ncx(''ncx_update_filevars(gcbo)'')');
  theVars3  = uicontrol('style','popupmenu','units','normalized','position', [x+2*w+d y-d-h       2*h+w     h      ],'string','no var','tag','ncx_varnames3',...
                        'callback',['ncx(''ncx_dim_ctrls(gcbo,3)'')']);
  frame32   = uicontrol('style','frame',    'units','normalized','position', [x+W-d/5       y-H     Wmax-W+d/5 H],'tag','ncx_frame32');

  is3pcolor = uicontrol('style','radiobutton','units','normalized','position', [x+W+d       y-h-d     2*w h],'string','pcolor','tag','ncx_3rdispcolor');
  is3cont   = uicontrol('style','radiobutton','units','normalized','position', [x+W+d       y-h-d-h   2*w h],'string','contour','tag','ncx_3rdiscont');
  set(is3pcolor,'callback','ncx(''rbuttons(gcbo,{''''ncx_3rdispcolor'''',''''ncx_3rdiscont''''})'')','value',0);
  set(is3cont,  'callback','ncx(''rbuttons(gcbo,{''''ncx_3rdispcolor'''',''''ncx_3rdiscont''''})'')','value',1);

  contVals3 = uicontrol('style','edit',       'units','normalized','position', [x+W+d       y-h-d-h-h 2*w h],'string','contours','tag','ncx_3rdcvals');
  contcolor3= uicontrol(                      'units','normalized','position', [x+W+d       y-4*h-d   0.5*w h],'string','c','callback','ncx(''color_button(gcbo)'')','tag','ncx_3rdccolor');
  contclab3 = uicontrol(                      'units','normalized','position', [x+W+d+w/2   y-4*h-d   1.5*w h],'string','clab','callback','ncx(''ncx_contlabel(gcbo,3)'')');

  frame3Handles = [frame3;theFiles3;theVars3;frame32;is3pcolor;is3cont;contVals3;contcolor3;contclab3];

% --------------------------------------------------------------------
% frame 4:
% --------------------------------------------------------------------
  y=y-H - d;
  frame4    = uicontrol('style','frame',    'units','normalized','position', [x       y-H     W H],'tag','ncx_frame4');
  theFiles4 = uicontrol('style','popupmenu','units','normalized','position', [x+d     y-d-h       2*w     h      ],'string','[0]','tag','ncx_files4',...
                        'callback','ncx(''ncx_update_filevars(gcbo)'')');
  theVars4  = uicontrol('style','popupmenu','units','normalized','position', [x+2*w+d y-d-h       2*h+w     h      ],'string','no var','tag','ncx_varnames4',...
    'callback',['ncx(''ncx_dim_ctrls(gcbo,4)'')']);
  txtX      = uicontrol('style','text',    'units','normalized','position', [x+2*w+d+2*h+w y-d-h       h     h      ],'string','X','FontWeight','bold');

  frame4Handles = [frame4;theFiles4;theVars4;txtX];

% --------------------------------------------------------------------
% frame 5:
% --------------------------------------------------------------------
  %y=y-H - d;
  x=x+W+d;
  frame5    = uicontrol('style','frame',    'units','normalized','position', [x       y-H     W H],'tag','ncx_frame5');
  theFiles5 = uicontrol('style','popupmenu','units','normalized','position', [x+d     y-d-h       2*w     h      ],'string','[0]','tag','ncx_files5',...
                        'callback','ncx(''ncx_update_filevars(gcbo)'')');
  theVars5  = uicontrol('style','popupmenu','units','normalized','position', [x+2*w+d y-d-h       2*h+w     h      ],'string','no var','tag','ncx_varnames5',...
    'callback',['ncx(''ncx_dim_ctrls(gcbo,5)'')']);
  txtY      = uicontrol('style','text',    'units','normalized','position', [x+2*w+d+2*h+w y-d-h       h     h      ],'string','Y','FontWeight','bold');

  frame5Handles = [frame5;theFiles5;theVars5;txtY];

% --------------------------------------------------------------------
% m_map:
% --------------------------------------------------------------------
if 0
  mmapCb   = uicontrol('style','checkbox','units','normalized','position', [x0 d 2.5*w h],'string','use m_map','tag','ncx_useM_map','value',0);
  strProj = {'Albers Equal-Area Conic','Lambert Conformal Conic','Stereographic','Orthographic','Azimuthal Equal-area','Azimuthal Equidistant','Gnomonic','Satellite','Mercator','Miller Cylindrical','Equidistant Cylindrical','UTM','Transverse Mercator','Sinusoidal','Gall-Peters','Hammer-Aitoff','Mollweide','Oblique Mercator'};
  mmapProj = uicontrol('style','popupmenu','units','normalized','position', [x0+2.5*w d 3*w h],'string',strProj,'tag','ncx_projM_map');
 mmapOptions = uicontrol(                  'units','normalized','position', [x0+5.5*w d 2*w h],'string','options','tag','ncx_optsM_map','callback','options_mmap(gcbo)');
end

% --------------------------------------------------------------------
% axes:
% --------------------------------------------------------------------
  xAX = x0+Wmax+5*d;
  yAX = 0.35;
  wAX = 1-xAX-1*d;
  hAX = 1-yAX-d;
  ax = axes('units','normalized','position',[xAX yAX wAX hAX]);

  wcb = wAX*.7;
  xcb = 1-d-wcb;
  ycb = 0.27;
  cbax = axes('units','normalized','position',[xcb ycb wcb h],'tag','ncx_colorbar'); axis off
  axes(ax);

% -------------------------etc

% zoom:
  y = 0.2;
  zoomInc = uicontrol(               'string','+','callback','ncx(''ncx_zoom(gcbo)'')','units','normalized','position',[1-d-h y     h h]);
  zoomVal = uicontrol('style','edit','string','0','callback','ncx(''ncx_zoom(gcbo)'')','units','normalized','position',[1-d-h y-h   h h],'tag','ncx_zoomVal');
  zoomDec = uicontrol(               'string','-','callback','ncx(''ncx_zoom(gcbo)'')','units','normalized','position',[1-d-h y-2*h h h]);
  zoomC   = uicontrol('style','checkbox','string','ref','callback','ncx(''ncx_zoom(gcbo)'')','units','normalized','position',[1-d-2*h y-3*h   2*h h],'tag','ncx_zoomCentre');

% caxis:
  theCaxis = uicontrol('style','edit','units','normalized','position', [1-d-2*h-2*w y 2*w h],'string','Caxis_','tag','ncx_Caxis');


% store frame's handles:
  setappdata(gcf,'frame_2',frame2Handles);
  setappdata(gcf,'frame_3',frame3Handles);
  setappdata(gcf,'frame_4',frame4Handles);
  setappdata(gcf,'frame_5',frame5Handles);


% align left popupmenu's strings
  pum=findobj(gcf,'style','popupmenu');
  set(pum,'HorizontalAlignment','left');

% ncx_mode:
  ncx_modes([],2);

% theme:
  ncx_theme

else
  set(gcf,'Pointer','arrow');
  eval(varargin{1})
end

function ncx_paths
base = fileparts(which(mfilename));
sep  = filesep;
pathData  = [base,sep,'data'];
pathUtils = [base,sep,'utils'];

% about netcdf
% matlab version: some paths are different in matlab < 6
v=version;
ir=find(v=='R');
ie=find(v==')');
v=str2num(v(ir+1:ie-1));
if isequal(computer,'PCWIN')
  dirs_netcdf = {
    [base,sep,'netcdf',sep,'netcdf'],
    [base,sep,'netcdf',sep,'netcdf',sep,'ncsource'],
    [base,sep,'netcdf',sep,'netcdf',sep,'nctype'],
    [base,sep,'netcdf',sep,'netcdf',sep,'ncutility'],
    [base,sep,'netcdf',sep,'mexcdf',sep,'pcwin'],
  };
else
  if v >=12
    dirs_netcdf = {
      [base,sep,'netcdf',sep,'netcdf'],
      [base,sep,'netcdf',sep,'netcdf',sep,'ncsource'],
      [base,sep,'netcdf',sep,'netcdf',sep,'nctype'],
      [base,sep,'netcdf',sep,'netcdf',sep,'ncutility'],
      [base,sep,'netcdf',sep,'mexcdf',sep,'v6_glx'],
    };
  else
    dirs_netcdf = {
      [base,sep,'netcdf',sep,'netcdf'],
      [base,sep,'netcdf',sep,'netcdf',sep,'ncsource'],
      [base,sep,'netcdf',sep,'netcdf',sep,'nctype'],
      [base,sep,'netcdf',sep,'netcdf',sep,'ncutility'],
      [base,sep,'netcdf',sep,'mexcdf',sep,'v5_glx'],
    };
  end
end

% add to path:
addpath(pathData);
addpath(pathUtils);
addNetCDF = 0;
if isempty(which('mexcdf')) | isempty(which('netcdf')) | isempty(which('mexcdf53'))
  addNetCDF = 1;
  for i=1:length(dirs_netcdf)
    addpath(dirs_netcdf{i});
  end

  % check if need to copy the dlls 
  if isequal(computer,'PCWIN') & v < 14 % not needed for R14!!
    if v >=12
      dllFolder = [matlabroot,filesep,'bin',filesep,'win32'];
    else
      dllFolder = [matlabroot,filesep,'bin'];
    end
    dllFiles = {'netcdf.dll','mexcdf53.dll'};

    for i=1:length(dllFiles)
      theSource = [dirs_netcdf{end},filesep,dllFiles{i}];
      if ~exist([dllFolder,filesep,dllFiles{i}])
        try
          copyfile(theSource,dllFolder);
        catch
          disp([':: ',mfilename,' : unable to copy ',theSource,' to ',dllFolder]);
        end
      end
    end

  end
end

% store the added path:
if addNetCDF
  addedPath = dirs_netcdf;
else
  addedPath = {};
end
addedPath{end+1} = pathData;
addedPath{end+1} = pathUtils;
setappdata(gcf,'addedPath',addedPath);

function ncx_close
addedPath = getappdata(gcf,'addedPath');
for i=1:length(addedPath)
  rmpath(addedPath{i});
end
closereq
