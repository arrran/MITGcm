function ncx_zoom(theHandle)

if nargin  == 0
  theHandle = findobj(gcf,'tag','ncx_zoomVal');
end
str = get(theHandle,'string');

if isequal(str,'ref')
  % show zoom centre:
  zoomC = findobj(gcf,'tag','ncx_zoomCentre');
  if get(zoomC,'value')
    zoomc('add')
  else
    cHandle = findobj(gcf,'tag','ncx_zoomC12');
    try
      delete(cHandle)
    end
  end

else
  % apply zoom

  valHandle = findobj(gcf,'tag','ncx_zoomVal');
  val = get(valHandle,'str'); val = str2num(val);

  %%  [rx1,rx2,ry1,ry2,rz1,rz2] = range_on_axis(gca,'Clipping','on');
  % use only region border:
  [rx1,rx2,ry1,ry2,rz1,rz2] = range_on_axis(gca,'tag','ncx_DataBorder');
  xlim([rx1 rx2]);
  ylim([ry1 ry2]);

  dx = diff(xlim)/2;
  dy = diff(ylim)/2;
  cHandle = findobj(gcf,'tag','ncx_zoomC12');
  try
    cHandle=cHandle(1);
    xc = get(cHandle,'xdata');
    yc = get(cHandle,'ydata');
  catch
    xc = sum(xlim)/2;
    yc = sum(ylim)/2;
  end

  switch str
    case {'+'}
      zoomVal = 1;
    case {'-'}
      zoomVal = -1;
    otherwise
      zoomVal = 0;
  end
  n=val+zoomVal;
  set(valHandle,'string',n);

  xlim([xc-dx/2^n xc+dx/2^n]);
  ylim([yc-dy/2^n yc+dy/2^n]);
end
