
clear;close;

[amp, header ]= sac2mat('E.sac');


dt=header(1);
npts=header(80);
b=header(6);
Fs=1/dt;
t= b : dt : (npts-1).*dt;
%figure(1)
%plot(t,amp)

% convert cm -> m
xdft = fft(amp*0.01);
xdft = xdft(1:npts/2+1);
psdx = (1/(Fs*npts)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 1./(0:Fs/length(amp):Fs/2);


% NLNM
[NLNM_data]=load('NLNM.txt');
period2=NLNM_data(:,1);
% convert to velocity
NLNM=NLNM_data(:,2)+NLNM_data(:,3).*log10(period2)+20*log10(period2/2/pi);

% NHNM
[NHNM_data]=load('NHNM.txt');
period1=NHNM_data(:,1);
% convert to velocity
NHNM=NHNM_data(:,2)+NHNM_data(:,3).*log10(period1)+20*log10(period1/2/pi);


semilogx(freq,10*log10(psdx))
hold on
plot(period2,NLNM)
plot(period1,NHNM)
grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
ylim([-220 -90])
xlim([0.01 10])
%set(gca, 'xdir','reverse')