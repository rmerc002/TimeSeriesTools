subLen_full = 300; 
maxWarp_full = 60;
numRandIters = 1000;
downsampleRatio = 10;

motionColumn = 64;
dataA_full = S1Drill{:,motionColumn};
dataB_full = S4Drill{:,motionColumn};

activityLabel = 407521;
maskA_full = S1Drill{:,250} == activityLabel;
maskB_full = S4Drill{:,250} == activityLabel;
classificationWindowA_full = zeros(size(maskA_full));

signalCount = 0;
startIndex = 1;
for i = 2:length(dataA_full)
    if maskA_full(i) > maskA_full(i-1)
        startIndex = i;
    elseif maskA_full(i) < maskA_full(i-1)
        if mod(signalCount,2) == 1
            maskA_full(startIndex:i) = 0;
        else
            classificationWindowA_full(startIndex-ceil(subLen_full/2):startIndex+ceil(subLen_full/2)) = 1;
        end
        signalCount = signalCount + 1;
    end
end



signalCount = 0;
startIndex = 1;
brange = nanstd(dataB_full);
for i = 2:length(dataB_full)
    if maskB_full(i) > maskB_full(i-1)
        startIndex = i;
    elseif maskB_full(i) < maskB_full(i-1)
        if mod(signalCount,2) == 1
            dataB_full(startIndex:i) = dataB_full(startIndex:i).*randn(i-startIndex+1,1);
        end
        signalCount = signalCount + 1;
    end
end

dataA = dataA_full(1:downsampleRatio:end);
dataB = dataB_full(1:downsampleRatio:end);

maskA = maskA_full(1:downsampleRatio:end);
maskB = maskB_full(1:downsampleRatio:end);
classificationWindowA = classificationWindowA_full(1:downsampleRatio:end);
groundTruthA = maskA;

subLenSeries = getSubLenSeries(40, 2000, 70);
indicatorHeatmap = nan(length(subLenSeries), length(dataA));
indices = ones(length(subLenSeries), length(dataA));
for subLenIndex = 1:length(subLenSeries)
    subLen_full = subLenSeries(subLenIndex);
    subLen = ceil(subLen_full/downsampleRatio)
    maxWarp = ceil(subLen*0.2);

    randRatio = 1/10;

    tic;
    [indicator] = chimp(dataA, dataB, subLen, maxWarp, 20);
    toc;
    indicatorHeatmap(subLenIndex,1:length(indicator)) = indicator;
    % figure; plot(distancesAB);
end