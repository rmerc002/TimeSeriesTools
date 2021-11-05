function [timeSeriesCombined, groundTruthCombined] = plotWeaklyLabeled(signals1,signals2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subLength = size(signals1,2);

noiseLength = 3*subLength;
[timeSeries1, groundTruth1] = generateSignalInRandomWalk(signals1,noiseLength);
[timeSeries2, groundTruth2] = generateSignalInRandomWalk(signals2,noiseLength);

timeSeriesCombined = [timeSeries1,timeSeries1(end)+timeSeries2];
groundTruthCombined = [groundTruth1,2*groundTruth2];

mp1 = mpx(timeSeries1',ceil(subLength/2),subLength);
mp1 = abs(mp1)';
mp2 = mpx(timeSeries2',ceil(subLength/2),subLength);
mp2 = abs(mp2)';
mpC = mpx(timeSeriesCombined',ceil(subLength/2),subLength);
mpC = abs(mpC)';

timeIndices = 1:length(timeSeriesCombined);


figure;
tiledlayout(2,1);
ax1 = nexttile;
hold on;
plot(timeSeriesCombined, 'Color', [0.5,0.5,0.5]);
coloredSignals = timeSeries1;
coloredSignals(groundTruth1~=1) = nan;
plot(coloredSignals,'g');
coloredSignals = timeSeries2;
coloredSignals(groundTruth2~=1) = nan;
plot((length(timeSeries1)+1:length(timeSeries1)+length(timeSeries2)),timeSeries1(end)+coloredSignals,'r');

for index=2:length(timeSeriesCombined)
   if groundTruthCombined(index)> groundTruthCombined(index-1)
       plot([index,index],[min(timeSeriesCombined),max(timeSeriesCombined)], 'Color', [0.9,0.9,0.9]);
   end
end

xlim([1,length(timeSeriesCombined)]);
hold off;

ax2 = nexttile;
hold on;
plot(mp1);
plot([nan(1,length(timeSeries1)),mp2]);
plot(mpC);

for index=2:length(timeSeriesCombined)
   if groundTruthCombined(index)> groundTruthCombined(index-1)
       plot([index,index],[min(mpC),max(mpC)], 'Color', [0.9,0.9,0.9]);
   end
end

ylim([-0.1,1.1]);
xlim([1,length(timeSeriesCombined)]);
hold off;

linkaxes([ax1 ax2], 'x');
end

