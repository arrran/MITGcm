function color_button(theHandle)


c=get(theHandle,'ForegroundColor');
c = uisetcolor(c);
set(theHandle,'ForegroundColor',c);
