function res=explode(str,s,n)
%EXPLODE   Split a string by string
%
%   Syntax:
%      OUT = EXPLODE(STR,S)
%
%   Inputs:
%      STR   The string
%      S     The needle
%      N     Substring indice to output   [ <none> ]
%
%   Output:
%      RES   Cell array if N is not specified; a string, otherwise
%
%   Examples:
%      str   = '09-Feb-2004';
%      year  = explode(str,'-',3);  % '2004'
%      parts = explode(str,'-');    % '09'    'Feb'    '2004'
%
%   MMA 9-6-2004, martinho@fis.ua.pt

%   Department of Physics
%   University of Aveiro, Portugal

%   **-08-2005 - Allow needle length > 1

if nargin == 3
   index=n;
else
   index=[];
end

i=findstr(str,s);
if isempty(i)
  res={str};
  return

elseif length(i) >=1
  res{1} = str(1:i-1);
  cont=1;
  if length(i) > 1
    for n=1:length(i)-1
      i1=i(n)+length(s);
      i2=i(n+1)-1;
      cont=cont+1;
      res{cont}=str(i1:i2);
    end
  end
  cont=cont+1;
  evalc('res{cont}=str(i(end)+length(s):end);','');
end

if ~isempty(index)
   res=res{index};
end
