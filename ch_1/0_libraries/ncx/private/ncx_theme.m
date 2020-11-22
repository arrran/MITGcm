function ncx_theme(theType)
%NCX_THEME   Apply a theme to the ncx gui
%
%   Syntax:
%      NCX_THEME(THETYPE)
%
%   Input:
%      THETYPE   'all' (default)
%                'figure'
%                'axes'
%                'edit'
%                'frame'
%                'checkbox'
%                'radiobutton'
%                'checkbox'
%                'listbox'
%                'text'
%                'pushbutton'
%                'togglebutton'
%                'popupmenu','
%                'chind'  % change indices, string='>','<'
%                'colorb' % color buttons string='c'

if nargin == 0
  theType = 'all';
end

% figure backgound:
FigureBg  = [171 176 205]/255;
% axes:
AxesBg    = [237 238 243]/255;
AxesFg    = [ 52  57  88]/255;
% frames:
FramesBg  = [206 209 225]/255;
FramesFg  = AxesFg;
% editable areas
EditBg    = [231 233 241]/255;
EditFg    = [111   0   0]/255;
% checkboxes:
CheckbBg  = FramesBg;
CheckbFg  = AxesFg;
% radiobuttons:
RadiobBg  = FramesBg;
RadiobFg  = AxesFg;
% listboxes:
ListbBg   = FramesBg;
ListbFg   = AxesFg;
% text areas:
TextaBg   = FramesBg;
TextaFg   = AxesFg;
% pushbuttons:
PushbBg   = FigureBg;
PushbFg   = [ 71  79 120]/255;
% tooglebuttons:
TooglebBg = PushbBg;
TooglebFg = PushbFg;
% popupmenus:
PopupmBg  = PushbBg;
PopupmFg  = PushbFg;
% special:
IncindFg  = [176   0   0]/255; % to use in >  and <
ColorbFg  = [  0   0   0];     % color buttons

% --------------------------------------------------------------------
% apply
% --------------------------------------------------------------------

if ismember(theType,{'figure' ,'all'})
    set(gcf,'color',FigureBg);
end
if ismember(theType,{'axes','all'})
    ax=findobj(gcf,'Type','axes');
    set(ax,'color',AxesBg,'xcolor',AxesFg,'ycolor',AxesFg,'zcolor',AxesFg);
end
if ismember(theType,{'frame','all'})
    obj = findobj(gcf,'style','frame');
    set(obj,'backgroundcolor',FramesBg,'foregroundcolor',FramesFg);
end
if ismember(theType,{'edit','all'})
    obj = findobj(gcf,'style','edit');
    set(obj,'backgroundcolor',EditBg,'foregroundcolor',EditFg);
end
if ismember(theType,{'checkbox','all'})
    obj = findobj(gcf,'style','checkbox');
    set(obj,'backgroundcolor',CheckbBg,'foregroundcolor',CheckbFg);
end
if ismember(theType,{'radiobutton','all'})
    obj = findobj(gcf,'style','radiobutton');
    set(obj,'backgroundcolor',RadiobBg,'foregroundcolor',RadiobFg);
end
if ismember(theType,{'listbox','all'})
    obj = findobj(gcf,'style','listbox');
    set(obj,'backgroundcolor',ListbBg,'foregroundcolor',ListbFg);
end
if ismember(theType,{'text','all'})
    obj = findobj(gcf,'style','text');
    set(obj,'backgroundcolor',TextaBg,'foregroundcolor',TextaFg);
end
if ismember(theType,{'pushbutton','all'})
    obj = findobj(gcf,'style','pushbutton');
    set(obj,'backgroundcolor',PushbBg);
    % here I must be caareful not to change fgcolor of color buttons:
    colorb = findobj(gcf,'style','pushbutton','string','c');
    obj = setdiff(obj,colorb);
    set(obj,'foregroundcolor',PushbFg);
end
if ismember(theType,{'togglebutton','all'})
    obj = findobj(gcf,'style','togglebutton');
    set(obj,'backgroundcolor',TooglebBg,'foregroundcolor',TooglebFg);
end
if ismember(theType,{'popupmenu','all'})
    obj = findobj(gcf,'style','popupmenu');
    set(obj,'backgroundcolor',PopupmBg,'foregroundcolor',PopupmFg);
end
if ismember(theType,{'chind','all'})
    obj1 = findobj(gcf,'string','>');
    obj2 = findobj(gcf,'string','<');
    obj=[obj1; obj2];
    set(obj,'foregroundcolor',IncindFg);
end
if ismember(theType,{'colorb','all'})
    obj = findobj(gcf,'style','pushbutton','string','c');
    set(obj,'foregroundcolor',ColorbFg);
end
