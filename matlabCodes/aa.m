function G = aa(ff)
[z,p,k] = ellip(8,0.001,80,2*pi*7570,'s');
G = zpktf(z,p,k*10^(0.001/20),ff) .* zpktf([],-2*pi*32768,2*pi*32768,ff);
end
