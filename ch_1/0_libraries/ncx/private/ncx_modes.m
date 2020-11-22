function ncx_modes(theHandle,theMode)

objs=findobj('tag','ncx_menu_modes');
if ~isempty(theHandle)
  set(objs,'checked','off');
  set(theHandle,'checked','on');
end

frame2Handles = getappdata(gcf,'frame_2');
frame3Handles = getappdata(gcf,'frame_3');
frame4Handles = getappdata(gcf,'frame_4');
frame5Handles = getappdata(gcf,'frame_5');

dimName2 = findobj(gcf,'tag','dim_name2');
dimCtrl2 = findobj(gcf,'tag','dim_ctrl2'); dim2Handles = [dimName2; dimCtrl2];

dimName3 = findobj(gcf,'tag','dim_name3');
dimCtrl3 = findobj(gcf,'tag','dim_ctrl3'); dim3Handles = [dimName3; dimCtrl3];

dimName4 = findobj(gcf,'tag','dim_name4');
dimCtrl4 = findobj(gcf,'tag','dim_ctrl4'); dim4Handles = [dimName4; dimCtrl4];

dimName5 = findobj(gcf,'tag','dim_name5');
dimCtrl5 = findobj(gcf,'tag','dim_ctrl5'); dim5Handles = [dimName5; dimCtrl5];


switch theMode
  case 1
    set(frame2Handles,'visible','off');
    set(frame3Handles,'visible','off');
    set(frame4Handles,'visible','off');
    set(frame5Handles,'visible','off');
    try, set(dim2Handles,'visible','off'); end
    try, set(dim3Handles,'visible','off'); end
    try, set(dim4Handles,'visible','off'); end
    try, set(dim5Handles,'visible','off'); end
  case 2
    set(frame2Handles,'visible','off');
    set(frame3Handles,'visible','off');
    set(frame4Handles,'visible','on');
    set(frame5Handles,'visible','on');
    try, set(dim2Handles,'visible','off'); end
    try, set(dim3Handles,'visible','off'); end
    try, set(dim4Handles,'visible','on'); end
    try, set(dim5Handles,'visible','on'); end
  case 3
    set(frame2Handles,'visible','on');
    set(frame3Handles,'visible','off');
    set(frame4Handles,'visible','on');
    set(frame5Handles,'visible','on');
    try, set(dim2Handles,'visible','on'); end
    try, set(dim3Handles,'visible','off'); end
    try, set(dim4Handles,'visible','on'); end
    try, set(dim5Handles,'visible','on'); end
  case 4
    set(frame2Handles,'visible','on');
    set(frame3Handles,'visible','on');
    set(frame4Handles,'visible','on');
    set(frame5Handles,'visible','on');
    try, set(dim2Handles,'visible','on'); end
    try, set(dim3Handles,'visible','on'); end
    try, set(dim4Handles,'visible','on'); end
    try, set(dim5Handles,'visible','on'); end
end
