function rbuttons(theHandle,tagList)

for i=1:length(tagList)
  rbHandle = findobj(gcf,'tag',tagList{i});
  set(rbHandle,'value',0);
end

set(theHandle,'value',1);
