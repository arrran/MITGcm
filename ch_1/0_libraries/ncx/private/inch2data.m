function theScale=inch2data(theScale)

OrigAxUnits = get(gca,'Units');
if OrigAxUnits(1:3) == 'nor'
   OrigPaUnits = get(gcf, 'paperunits');
   set(gcf, 'paperunits', 'inches');
   figposInches = get(gcf, 'paperposition');
   set(gcf, 'paperunits', OrigPaUnits);
   axposNor = get(gca, 'position');
   axWidLenInches = axposNor(3:4) .* figposInches(3:4);
else
   set(gca, 'units', 'inches');
   axposInches = get(gca, 'position');
   set(gca, 'units', OrigAxUnits);
   axWidLenInches = axposInches(3:4);
end

scX = diff(get(gca, 'XLim'))/axWidLenInches(1);
scY = diff(get(gca, 'YLim'))/axWidLenInches(2);
sc = max([scX;scY]);


theScale = sc/theScale;
