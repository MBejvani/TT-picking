function [out, index] = zcr(in, slope, level)
% [out, index] = zcr(in, slope, level)
%
% Zero crossing rate of input data (x-axis crossings)
%
% out: Number of zero crossings
% index: indizes into matrix where zero-xings are found
% in:  data
% slope: (string, only first char is taken into account) 
%        'positive': all zero-xings with rising slopes are counted
%        'negative': all zero-xings with falling slopes
%        'all'     : all zero-xings (= default)
% level: (scalar) crossings over y = level (default: level = 0, i.e. zero-crossings)
%        
% Armin Günter
% Rev.1.0, 6.10.97
% Rev.2.0, 19.12.97: Input Argument check, Matrix support
% Rev.3.0, 6.3.98: optional index output, optional slope selection,
%    ending pts are not counted any more!

% Fehler abfangen:
error(nargchk(1,3,nargin))
if (isempty(in))
   y = [];
   return
end
error(nargchk(0,2,nargout))
if ~exist('level')
   level = 0;
end

temp = sign(in(:) - level(1));
% Eliminate zeros; otherwise a -x,0,y axis crossing would be counted twice! 
temp(temp==0) = 1;
if exist('slope')
   switch slope(1)
   case 'p'
      index = find(diff(temp) > 0);
   case 'n'
      index = find(diff(temp) < 0);
   otherwise
      index = find(diff(temp));
   end
else
   index = find(diff(temp));
end
out = length(index);