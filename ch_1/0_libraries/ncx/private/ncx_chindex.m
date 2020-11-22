function ncx_chindex(theHandle,theAction)

otherHandles = get(theHandle,'userdata');
dim  = otherHandles(1); dimn  = get(dim,'string');
ind1 = otherHandles(2); nind1 = str2num(get(ind1,'string'));
ind2 = otherHandles(3); nind2 = str2num(get(ind2,'string'));
ind3 = otherHandles(4); nind3 = str2num(get(ind3,'string'));
imax = otherHandles(5); nimax = str2num(get(imax,'string'));

switch theAction
  case {'inc'}
    if nind1==nind3
      nind1=nind1+nind2;
      if nind1>nimax
        nind1=1;
      end
    end
  case {'dec'}
    if nind1==nind3
      nind1=nind1-nind2;
      if nind1<1
        nind1=nimax;
      end
    end
  case {'manual'}
    % here I could check if the inserted values are correct... TODO
end

if ismember(theAction,{'inc','dec'})
  nind3=nind1;
  set(ind1,'string',nind1);
  set(ind3,'string',nind3);
end

% --------------------------------------
% change also other indices!! 2nd, 3rd var, X and Y
for i=2:5
  dim_nameTag      = ['dim_name',num2str(i)];
  objs = findobj(gcf,'tag',dim_nameTag);
  if ~isempty(objs)
    dims = get(objs,'string');
    % find the dim number (position)
    idim=[];
    if iscell(dims)
      for nd = 1:length(dims)
        if isequal(dims{nd},dimn)
          idim=nd;
          break
        end
      end
    else
      if isequal(dims,dimn)
        idim=1;
      end
    end
    if ~isempty(idim)
      % the dim was found, so change it:
      otherHandles = get(objs(idim),'userdata');

      ind1 = otherHandles(1); set(ind1,'string',nind1);
      ind2 = otherHandles(2); set(ind2,'string',nind2);
      ind3 = otherHandles(3); set(ind3,'string',nind3);
    end
  end
end

% disp if check box on:
cb = findobj(gcf,'tag','ncx_dispcb');
if get(cb,'value') & ~isequal(theAction,'manual')
  ncx('ncx_plot');
end
