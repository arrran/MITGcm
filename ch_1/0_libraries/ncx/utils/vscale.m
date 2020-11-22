function [val,n,ex] = vscale(M)
%VSCALE   Base 10 scale - first non zero algarism
%   Returns  the first non zero algarism of M or mean(mean(M)).
%
%   Syntax:
%      [VAL,IVAL,EX] = VSCALE(M)
%
%   Input:
%      M   Array 2D or less
%
%   Output:
%      VAL   Mean M with only one non zero algarism
%      IVAL  The non zero algarism
%      EX    The base 10 power
%
%   Example2:
%     [v,iv,e] = vscale(pi);       % returns   3, 3,  0
%     [v,iv,e] = vscale(pi*100);   % returns 300, 3,  2
%     [v,iv,e] = vscale(rand(10)); % returns 0.5, 5, -1
%
%   MMA 6-8-2005, martinho@fis.ua.pt

%   Department of Physics
%   University of Aveiro, Portugal

if prod(size(M)) ~= 1
  M=mean(mean(M));
end

str = sprintf('%e',M);
tmp = explode(str,'e');
n   = round(str2num(tmp{1}));
ex  = str2num(tmp{2});
val = n*10.^ex;
