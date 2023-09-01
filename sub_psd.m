function [ psdx ,freq ] = sub_psd(amp,dt,npts)

Fs=1/dt;
xdft = fft(amp);
xdft = xdft(1:npts/2+1);
psdx = (1/(Fs.*npts)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 1./(0:Fs/length(amp):Fs/2);
psdx=10*log10(psdx);

end
%set(gca, 'xdir','reverse')