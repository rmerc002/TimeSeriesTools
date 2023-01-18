mpAB = mpx_ABBA_v2(ts, tsB, mm)/sqrt(4*mm);

minLength = min(size(frequencyProfile,2), length(mpAB));
frequencyCP = frequencyProfile(1:minLength);

for ii = 1:numThresholds
    frequencyCP(ii,:) = frequencyProfile(ii,1:minLength).*mpAB(1:minLength)';
%     frequencyCP(ii,:) = frequencyCP(ii,:)./max(frequencyCP(ii,:));
end


% figure;
% imagesc(frequencyCP);
% set(gca,'YDir','normal');


figure;
tiledlayout(2,1);
ax1 = nexttile();
plot(ts);

ax2 = nexttile();
imagesc(frequencyCP);
set(gca,'YDir','normal')
 colormap(gray);
 
linkaxes([ax1, ax2], 'x');
