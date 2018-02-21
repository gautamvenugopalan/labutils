function ftf = FotonFilter(fname,modname,fm_num,ff,varargin)
% FOTONFILTER takes a FOTON filter file, the module name and the filter
% number and returns the TF
%
% use 'DAC' as 5th argument to add the DAC sample-and-hold rolloff
%
% Ex 1:  darm_fm1 = FotonFilter('H2LSC.txt','DARM',1,f);
%

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

pp = readFilterFile(fname);



% S-domain stuff
  fs = eval(['pp.' modname '(' num2str(fm_num) ').fs']);
  digw = (2 * pi / fs) .* ff;    % Convert from Hz to rad/sample
  s = exp(-i * digw);           % S for this calculation
  ss = s .* s;                  % S^2



sc = eval(['pp.' modname '(' num2str(fm_num) ').soscoef']);

h = ones(size(ff));

for m = 1:size(sc, 1)
   num = sc(m, 1) + sc(m, 2) .* s + sc(m, 3) .* ss;
   den = sc(m, 4) + sc(m, 5) .* s + sc(m, 6) .* ss;
   h = h .* num ./ den;
end

switch upper(dacw)
 case 'DAC'
  d2a = (sin(pi*ff/fs)./(pi*ff/fs)).*exp(-i*pi*ff/fs);
  ftf = h.*d2a;
 otherwise
  ftf = h;
end

