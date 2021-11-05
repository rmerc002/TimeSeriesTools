meancup = nanmean(S2Drill{:,2});
stdcup = nanstd(S2Drill{:,2});
cup = (S2Drill{:,2}-meancup)/stdcup;
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
groundTruthIndicesTest = [];
i = 1;
signalLength = 150;
while i <= length(isDrinkingFilteredNeg)
    if isDrinkingFiltered(i) == true
        groundTruthIndicesTest(end+1) = i;
        i = i + ceil(signalLength*5);
    else
        i = i + 1;
    end
end

testTS = S2Drill{:,1};
testMean = nanmean(testTS);
testTS(isnan(testTS)) = testMean;

tiledlayout(5,1);
ax1 = nexttile;
plot(testTS);
title("Test Data");

ax2 = nexttile;
plot(cup);

ax3 = nexttile;
plot(isDrinkingRaw);

ax4 = nexttile;
plot(isDrinkingFiltered);

% ax5 = nexttile;
% plot(isDrinkingFilteredNeg);
% plot(isDrinkingFilteredPos);

ax5 = nexttile;
plot(0,0);
hold on;
for i = 1:length(groundTruthIndicesTest)
    gti = groundTruthIndicesTest(i);
    plot([gti,gti],[0,1],'g');
end
hold off;
title("Ground Truth Indices Test");

linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');


%%% Output variables
figureSavePath = "C:\Users\Ryan\Dropbox\Education\UCR\Research\0_Projects\ActiveResearch\shapeletDiscovery\Experiments\2020-11-16_Opportunity\Results";
[tpIndices, fpIndices, fnIndices] = evaluateAccuracy(testTS, signalLength, groundTruthIndicesTest, queryShapelet, queryThreshold, thresholdConfidence, figureSavePath);