function [timeSeriesAB, groundTruthAB] = plotWeaklyLabeled_V03(signalsA,signalsB, numberNearestNeighbors)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

subLength = size(signalsA,2);
% numberNearestNeighbors = 1;

noiseLength = 3*subLength;
[timeSeriesA, groundTruthA] = generateSignalInRandomWalk(signalsA,noiseLength);
[timeSeriesB, groundTruthB] = generateSignalInRandomWalk(signalsB,noiseLength);
% %temporary
% timeSeriesB(1:2:end) = 1;
% timeSeriesB(2:2:end) = 0;

timeSeriesAB = [timeSeriesA,timeSeriesA(end)+timeSeriesB];
groundTruthAB = [groundTruthA,2*groundTruthB];

mpAA = mpx_knn(timeSeriesA',ceil(subLength/2),subLength,numberNearestNeighbors);
mpAA = 1 - abs(mean(mpAA,1));
mpBB = mpx_knn(timeSeriesB',ceil(subLength/2),subLength,numberNearestNeighbors);
mpBB = 1 - abs(mean(mpBB,1));
mpAB = mpx_AB_knn(timeSeriesA',timeSeriesB',subLength,numberNearestNeighbors);
mpAB = 1 - abs(mean(mpAB,1));
mpBA = mpx_AB_knn(timeSeriesB',timeSeriesA',subLength,numberNearestNeighbors);
mpBA = 1 - abs(mean(mpBA,1));

timeIndices = 1:length(timeSeriesAB);


figure;
tiledlayout(2,1);
ax1 = nexttile;
hold on;
plot(timeSeriesAB, 'Color', [0.5,0.5,0.5]);
coloredSignals = timeSeriesA;
coloredSignals(groundTruthA~=1) = nan;
plot(coloredSignals, 'Color', [0,0.9,0]);
coloredSignals = timeSeriesB;
coloredSignals(groundTruthB~=1) = nan;
plot((length(timeSeriesA)+1:length(timeSeriesA)+length(timeSeriesB)),timeSeriesA(end)+coloredSignals,'Color',[0,0.8,1]);

for index=2:length(timeSeriesAB)
   if groundTruthAB(index)> groundTruthAB(index-1)
       plot([index,index],[min(timeSeriesAB),max(timeSeriesAB)], 'Color', [0.75,0.75,0.75]);
   end
end

xlim([1,length(timeSeriesAB)]);
hold off;

ax2 = nexttile;
hold on;
plot(mpAA,'Color',[0,0.9,0]);
plot(mpAB,'Color',[0,0.5,0]);
plot(max(0,mpAA-mpAB),'Color',[1,0,0]);
plot([nan(1,length(timeSeriesA)),mpBB],'Color',[0,0.8,0.8]);
plot([nan(1,length(timeSeriesA)),mpBA],'Color',[0,0.3,1]);
plot([nan(1,length(timeSeriesA)),max(0,mpBB-mpBA)],'Color',[1,0,0]);

for index=2:length(timeSeriesAB)
   if groundTruthAB(index)> groundTruthAB(index-1)
       plot([index,index],[0,1], 'Color', [0.75,0.75,0.75]);
   end
end

ylim([-0.1,1.1]);
xlim([1,length(timeSeriesAB)]);
hold off;

linkaxes([ax1 ax2], 'x');
end

