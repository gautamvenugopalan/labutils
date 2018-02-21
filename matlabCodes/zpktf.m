function ftf = zpktf(zero,peak,gain,ff)
    ftf = ones(size(ff));
    for n = 1 : length(zero)
    ftf = ftf.*(2 * pi * 1i * ff - ones(size(ff)) * zero(n));
    end
    for n = 1 : length(peak)
    ftf = ftf./(2 * pi * 1i * ff - ones(size(ff)) * peak(n));
    end
    ftf = gain * ftf;
end