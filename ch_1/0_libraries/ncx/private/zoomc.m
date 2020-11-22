function zoomc(what,pos)

cp = get(gca,'CurrentPoint');
x0 = cp(1,1);
y0 = cp(1,2);

if isequal(what,'add')
  if nargin == 2
    xc=pos(1);
    yc=pos(2);
  else
    xc = sum(xlim)/2;
    yc = sum(ylim)/2;
  end

  ish=ishold; hold on
  p1=plot(xc,yc,'+'); set(p1,'MarkerSize',20,'tag','ncx_zoomC12','erasemode','xor');
  p2=plot(xc,yc,'s'); set(p2,'MarkerSize',10,'tag','ncx_zoomC12','erasemode','xor');
  if ~ish, hold off; end
  set(gcf,'WindowButtonDownFcn','ncx(''zoomc(''''get'''')'')');
elseif isequal(what,'get')
  tag = get(gco,'tag');
  if isequal(tag,'ncx_zoomC12')
    set(gcf,'WindowButtonMotionFcn','ncx(''zoomc(''''move'''')'')');
    set(gcf,'WindowButtonUpFcn','ncx(''zoomc(''''stop'''')'')');
  end
elseif isequal(what,'move')
  cHandle = findobj(gcf,'tag','ncx_zoomC12');
  set(cHandle,'xdata',x0,'ydata',y0);
elseif isequal(what,'stop')
  % put cHandle at  the centre:
  cHandle = findobj(gcf,'tag','ncx_zoomC12'); cHandle=cHandle(1);
  xc = get(cHandle,'xdata');
  yc = get(cHandle,'ydata');
  dx = diff(xlim);
  dy = diff(ylim);
  xlim([xc-dx/2 xc+dx/2])
  ylim([yc-dy/2 yc+dy/2])

  set(gcf,'WindowButtonMotionFcn','');
  set(gcf,'WindowButtonUpFcn','');
end
