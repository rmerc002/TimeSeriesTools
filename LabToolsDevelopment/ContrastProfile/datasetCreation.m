filename = "";
type = "segmented";
needsTranspose = false;
classIndex = 1;

numDataSamples = 5;

data = sortrows(CBFTRAIN{:,:},classIndex);
classes = data(:,classIndex);
data(:,classIndex) = [];

figure;
numPlotSamples = 50;
numSamples = size(data,1);
step = floor(numSamples/numPlotSamples);
plot(0,0);
hold on;
for plotIndex = 1:numDataSamples
   dataIndex = plotIndex*step;
   tempTS = data;
   tempMin = min(tempTS);
   tempMax = max(tempTS);
   tempRange = max(1e-5, tempMax-tempMin);
   tempTSNorm = (tempTS-tempMin)/tempRange;
   
   plot(-plotIndex+0.9*tempTSNorm);
end
hold off;


positiveSamples = data(1:numDataSamples,:);
tsBaseLength = numDataSamples*2*10000;
randomTS = getRandWalk(tsBaseLength);

signalLength = size(positiveSamples,2);
subLength = signalLength;
sampleSectionRange = tsBaseLength/numDataSamples;
randOffset = randi(sampleSectionRange, numDataSamples,1)';
boundedRandIndices = randOffset + ((1:numDataSamples)-1)*sampleSectionRange;

positiveTS = [0];
for i = 1:numDataSamples
    ro = randOffset(i);
    part1 = zscore(getRandWalk(ro));
    part1 = part1 - part1(1) + positiveTS(end);
    part2 = zscore(positiveSamples(i,:));
    part2 = part2 - part2(1) + part1(end);
    part3 = zscore(getRandWalk(sampleSectionRange - ro));
    part3 = part3 - part3(1) + part2(end);
    positiveTS = [positiveTS, part1, part2 ,part3];
end

negativeTS = getRandWalk(length(positiveTS));

groundTruthIndices = 1 + boundedRandIndices + ((1:numDataSamples)-1)*signalLength;

figure;

plot(positiveTS);
hold on;
for i = 1:numDataSamples
    index = groundTruthIndices(i);
   plot(index:index+subLength-1, positiveTS(index:index+subLength-1),'green');
end
hold off;

