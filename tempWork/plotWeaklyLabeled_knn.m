function [timeSeriesCombined, groundTruthCombined] = plotWeaklyLabeled_knn(signals1,signals2, numberNearestNeighbors)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
subLength = size(signals1,2);
% numberNearestNeighbors = 1;
noiseLength = 3*subLength;
[timeSeries1, groundTruth1] = generateSignalInRandomWalk(signals1,noiseLength);
[timeSeries2, groundTruth2] = generateSignalInRandomWalk(signals2,noiseLength);
%temporary
timeSeries2(1:2:end) = 1;
timeSeries2(2:2:end) = 0;

timeSeriesCombined = [timeSeries1,timeSeries1(end)+timeSeries2];
groundTruthCombined = [groundTruth1,2*groundTruth2];

[mp,mpi] = mpx_knn(timeSeriesCombined, ceil(subLength/2),subLength,numberNearestNeighbors);
mp_locations = mpi;
mp_locations(:,1:end/2) = mp_locations(:,1:end/2) < length(timeSeriesCombined)/2;
mp_locations(:,end/2:end) = mp_locations(:,end/2:end) > length(timeSeriesCombined)/2;
mp_locations = mp_locations*2-1;
mp_locations = mp_locations.*(1-mp); %new line to weight matches
mpi_sums = sum(mp_locations,1);

timeIndices = 1:length(timeSeriesCombined);

figure;
tiledlayout(2,1);
ax1 = nexttile;
hold on;
plot(timeSeriesCombined, 'Color', [0.5,0.5,0.5]);
coloredSignals = timeSeries1;
coloredSignals(groundTruth1~=1) = nan;
plot(coloredSignals,'Color',[0,0.9,0]);
coloredSignals = timeSeries2;
coloredSignals(groundTruth2~=1) = nan;
plot((length(timeSeries1)+1:length(timeSeries1)+length(timeSeries2)),timeSeries1(end)+coloredSignals,'Color',[0,0.8,.8]);

for index=2:length(timeSeriesCombined)
   if groundTruthCombined(index)> groundTruthCombined(index-1)
       plot([index,index],[min(timeSeriesCombined),max(timeSeriesCombined)], 'Color', [0.9,0.9,0.9]);
   end
end
xlim([1,length(timeSeriesCombined)]);
hold off;

ax2 = nexttile;
hold on;
plot(mpi_sums,'r');
for index=2:length(timeSeriesCombined)
   if groundTruthCombined(index)> groundTruthCombined(index-1)
       plot([index,index],[-1*numberNearestNeighbors,numberNearestNeighbors], 'Color', [0.9,0.9,0.9]);
   end
end
xlim([1,length(timeSeriesCombined)]);
ylim([0,numberNearestNeighbors]);
hold off;

linkaxes([ax1 ax2], 'x');
end

