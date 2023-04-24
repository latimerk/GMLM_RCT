figure(1);
clf

pp = 1:57;

tts = NC.tts(NC.rr);
for ii = 1:5
    subplot(2,3,ii)
    hold on
    plot(tts, squeeze(nanmean(NC.V2( NC.rr,pp,ii),2)),'k')
    plot(tts, squeeze(nanmean(NC.W2( NC.rr,pp,ii),2)))
    plot(tts, squeeze(nanmean(LC.W2(NC.rr,pp,ii),2)))
    hold off
end