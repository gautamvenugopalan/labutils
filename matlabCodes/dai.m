

function G = dai(f)
G = sos2tf(readAIAAcoefs(),f,'DAC');
end

function soscoef = readAIAAcoefs(n)

if(nargin == 0)
    n = 2;
end
% open file
  fid = fopen('/opt/rtcds/rtscore/release/src/fe/controller.c');
  if( fid == -1 )
    errfmt('Unable to open filter file "%s".\n', fileName);
  end

  % init
  i = 1;
  % loop through file
  str = fgetl(fid);
  while( isstr(str) )

    % parse line into cell array of strings
    if( isempty(str) )
      arg = {};
    else
      arg = strread(str,'%s','delimiter',',');
    end

    % check line
    if( isempty(arg) )
      % empty or blank line
    elseif( strcmp(arg{1},'static double __attribute__ ((unused)) feCoeff4x[9] =' ) )
    %elseif( strcmp(arg{1}, 'static double feCoeff4x[9] =') )
        if(i == n)
        str = fgetl(fid);
        gain = strread(str,'%f','delimiter',['{',',']);
        str = fgetl(fid);
        coeff1 = strread(str,'%f','delimiter',',');
        str = fgetl(fid);
        coeff2 = strread(str,'%f','delimiter',['};',',']);
        i=i+1;
        gain = gain(2);
        coef = [coeff1 coeff2(1:4)]';
        soscoef = zeros(2, 6);
        for n = 1:2
          soscoef(n, :) = [1 coef(n, 3) coef(n,4) 1 coef(n,1) coef(n,2)];
        end
        soscoef(1, :) = soscoef(1, :) .* [gain gain gain 1 1 1];
        else
            i=i+1;
        end
    else
        
    end

    % read next line
    str = fgetl(fid);
  end

  fclose(fid);
  
  
end

function h = sos2tf(soscoefs,ff,dacw)
% S-domain stuff
  fs = 16384 *4;
  digw = (2 * pi / fs) .* ff;    % Convert from Hz to rad/sample
  s = exp(-1i * digw);           % S for this calculation
  ss = s .* s;                  % S^2


h = ones(size(ff));

for m = 1:size(soscoefs, 1)
   num = soscoefs(m, 1) + soscoefs(m, 2) .* s + soscoefs(m, 3) .* ss;
   den = soscoefs(m, 4) + soscoefs(m, 5) .* s + soscoefs(m, 6) .* ss;
   h = h .* num ./ den;
end

switch upper(dacw)
 case 'DAC'
  d2a = (sin(pi*ff/fs)./(pi*ff/fs)).*exp(-i*pi*ff/fs);
  ftf = h.*d2a;
 otherwise
  ftf = h;
end
end
