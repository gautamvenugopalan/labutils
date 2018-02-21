function G = ai(ff)
[z,p,k] = ellip(4,4,60,2*pi*7570,'s');
G = zpktf(z,p,k*10^(4/20),ff) .* zpktf([],-2*pi*13e3,2*pi*13e3,ff);
end
