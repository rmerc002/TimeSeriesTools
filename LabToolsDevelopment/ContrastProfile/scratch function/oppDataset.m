rng(1);
meancup = nanmean(S1Drill{:,2});
stdcup = nanstd(S1Drill{:,2});
cup = (S1Drill{:,2}-meancup)/stdcup;
isDrinkingRaw = abs(cup) > 1.5;

isDrinkingFiltered = false(length(isDrinkingRaw),1);
i = 1;
signalLength = 150;
while i <= length(isDrinkingFiltered)
    if isDrinkingRaw(i) == true
        startIndex = i;
        endIndex = min(length(isDrinkingFiltered), i + signalLength);
        isDrinkingFiltered(startIndex:endIndex) = true;
%         fprintf("startIndex: %d\n",startIndex);
        i = i + ceil(signalLength*1.5);
    else
        i = i + 1;
    end
end

isDrinkingFiltered(54230:end) = false;


isDrinkingFilteredNeg = isDrinkingFiltered;
isDrinkingFilteredPos = isDrinkingFiltered;
groundTruthIndicesTrain = [];
startEndIndicesNeg = [];
i = 1;
signalLength = 150;
while i <= length(isDrinkingFilteredNeg)
    if isDrinkingFiltered(i) == true
        startIndex = i;
        endIndex = min(length(isDrinkingFilteredNeg), i + signalLength);
        isDrinkingFilteredNeg(startIndex:endIndex) = false;
        isDrinkingFilteredPos(startIndex+200:endIndex+450) = false;
%         fprintf("startIndex: %d\n",startIndex);
        groundTruthIndicesTrain(end+1) = i;
        startEndIndicesNeg(end+1,:) = [startIndex, endIndex];
        i = i + ceil(signalLength*4);
    else
        i = i + 1;
    end
end


tiledlayout(4,1);
ax1 = nexttile;
plot(cup);

ax2 = nexttile;
plot(isDrinkingRaw);

ax3 = nexttile;
plot(isDrinkingFiltered);

ax4 = nexttile;
plot(isDrinkingFilteredNeg);
plot(isDrinkingFilteredPos);

linkaxes([ax1 ax2 ax3 ax4], 'x');

startIndexTrainNeg = 1;
endIndexTrainNeg = ceil(size(S1Drill,1)/2);

startIndexTrainPos = endIndexTrainNeg+1;
endIndexTrainPos = size(S1Drill,1);



trainNegative = zscore(S1Drill{startIndexTrainNeg:endIndexTrainNeg,1});
trainMean = nanmean(trainNegative);
trainSTD = nanstd(trainNegative);
trainNegative(isnan(trainNegative)) = trainMean;
for i = 1:size(startEndIndicesNeg,1)
    startIndex = startEndIndicesNeg(i,1);
    endIndex = startEndIndicesNeg(i,2);
    if startIndex <= endIndexTrainNeg  && endIndex <= endIndexTrainNeg
        fillerLength = endIndex - startIndex + 1;
        trainNegative(startIndex:endIndex) = (zscore(getRandWalk(fillerLength)) ) * trainSTD;
    end
end


trainPositive = S1Drill{startIndexTrainPos:endIndexTrainPos,1};
trainMean = nanmean(trainPositive);
trainPositive(isnan(trainPositive)) = trainMean;

groundTruthIndicesTrain = groundTruthIndicesTrain(groundTruthIndicesTrain >= startIndexTrainPos);
groundTruthIndicesTrain = groundTruthIndicesTrain - length(trainNegative);
figureSavePath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\shapeletDiscovery\Experiments\2020-11-16_Opportunity\Results";
[queryShapelet, queryThreshold, thresholdConfidence] = shapeletDiscovery(trainPositive, trainNegative, signalLength, groundTruthIndicesTrain, figureSavePath);