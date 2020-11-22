function ncx_coastl(theHandle)

coastlineButtonTag  = 'ncx_addCoastline';
coastlinePlottedTag = 'ncx_coastline';

theAction = '';
if isstr(theHandle)
  theAction = theHandle;
  theHandle = findobj(gcf,'tag',coastlineButtonTag);
end
addCL = get(theHandle,'checked');

if isequal(theAction,'add') & isequal(addCL,'on')
  cl = load('wcl_0_360.dat');
  % if x<0, then rot_longitude:
  % get it fro data border:
  dataBorder=findobj(gcf,'tag','ncx_DataBorder');
  try
    xborder = get(dataBorder,'xdata');
    xmin = min(xborder);
  catch
    xmin=0;
  end  
  clx = cl(:,1);
  cly = cl(:,2);
  if xmin<0
    [clx,cly] = rot_longitude(clx,cly);  
  end
  ish=ishold; hold on
  plot(clx,cly,'tag',coastlinePlottedTag);
  if ~ish, hold off; end

elseif isequal(theAction,'remove')
  clHandle = findobj(gca,'tag',coastlinePlottedTag);
  try, delete(clHandle); end

elseif isequal(theAction,'')
  if isequal(addCL,'on')
    set(theHandle,'checked','off');
    ncx_coastl('remove');
  else
    set(theHandle,'checked','on');
    ncx_coastl('add');
  end
end
