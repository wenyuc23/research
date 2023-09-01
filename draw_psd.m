function figure(1) = draw_psd(freq,psdave,period2,NLNM,period1,NHNM,tt)

figure (1); % PSD
semilogx(freq,psdave)
hold on
plot(period2,NLNM,'k')
plot(period1,NHNM,'k')
text(0.05,-91,'NHNM')
text(0.05,-167,'NLNM')
grid on
title(['PSD-' tt])
xlabel('Period (sec)')
ylabel('PSD (dB)')
ylim([-200 -80])
xlim([0.01 1])

end