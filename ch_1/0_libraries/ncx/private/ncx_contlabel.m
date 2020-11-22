function ncx_contlabel(theHandle,theVar)

switch theVar
  case 1
    vals = getappdata(gcf,'contours_1stvar');
  case 2
    vals = getappdata(gcf,'contours_2ndvar');
  case 3
    vals = getappdata(gcf,'contours_3rdvar');
end

cs=nan;
ch=nan;
try
  cs = vals.cs;
  ch = vals.ch;
end
if ishandle(ch)
  clabel(cs,'manual');
end

% in some cases view~=2 occurs after the clabel !!?, so:
view(2)
