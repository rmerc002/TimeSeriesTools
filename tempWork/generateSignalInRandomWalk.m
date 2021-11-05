function [synthTimeSeries,groundTruthLabels] = generateSignalInRandomWalk(signals,lengthWalk)
%assume signals is shape n,m (samples, subsequenceLength)
randSeed = mod(ceil(100000*sum(abs(signals(1,:))))-1,2^23);
rng(randSeed);
synthTimeSeries = zscore(getRandWalkBasic(lengthWalk));
groundTruthLabels = zeros(1,lengthWalk);
for sampleIndex = 1:size(signals,1)
    rng(mod(randSeed+sampleIndex,2^32));
    normalizedSignal = signals(sampleIndex,:);
    normalizedSignal = (normalizedSignal-normalizedSignal(1))/std(normalizedSignal);
    synthTimeSeries = [synthTimeSeries, (synthTimeSeries(end) + normalizedSignal)];
    tempRand = getRandWalkBasic(lengthWalk);
    tempRand = (tempRand - tempRand(1))/std(tempRand);
    synthTimeSeries = [synthTimeSeries, (synthTimeSeries(end) + tempRand)];
    
    groundTruthLabels = [groundTruthLabels, ones(1,size(signals,2)),zeros(1,lengthWalk)];
end
end

