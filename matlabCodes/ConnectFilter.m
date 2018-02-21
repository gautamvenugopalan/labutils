
function ftf = ConnectFilter(filename, modname,fm_numList,ff,varargin)
    
dacw = 0;
if nargin < 4
  error('Needs at least 4 arguments, see the help.')
  if nargin > 4
    dacw = varargin{1};
    if nargin > 5
      error('Too Many Arguments.')
    end
  end
end

pp = readFilterFile(filename);
fs = eval(['pp.' modname '(' num2str(fm_numList(1)) ').fs']);

h = ones(size(ff));

for n = 1:length(fm_numList)
    h = h .* FotonFilter(filename, modname,fm_numList(n),ff);
end


switch upper(dacw)
 case 'DAC'
  d2a = (sin(pi*ff/fs)./(pi*ff/fs)).*exp(-1i*pi*ff/fs);
  ftf = h.*d2a;
 otherwise
  ftf = h;
end
end
