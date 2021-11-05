subLen_full = 300; 
maxWarp_full = 60;
numRandIters = 1000;
downsampleRatio = 10;

% motionColumn = 231;
dataA_full = S1Drill{:,motionColumn};
dataA_full = dataA_full + randn(size(dataA_full))*1e-3;
dataB_full = S2Drill{:,motionColumn};
dataB_full = dataB_full + randn(size(dataB_full))*1e-3;
% dataA_full = getRandWalk(length(S1Drill{:,motionColumn}));
% dataB_full = getRandWalk(length(S2Drill{:,motionColumn}));

activityLabel = 407521;
maskA_full = S1Drill{:,250} == activityLabel;
maskB_full = S2Drill{:,250} == activityLabel;
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
%             classificationWindowA_full(startIndex-ceil(subLen_full/2):startIndex+ceil(subLen_full/2)) = 1;
            classificationWindowA_full(startIndex-ceil(subLen_full*2):startIndex+ceil(subLen_full*4)) = 1;

        end
        signalCount = signalCount + 1;
    end
end



signalCount = 0;
startIndex = 1;
brange = nanstd(dataB_full);
meanB = nanmean(dataB_full);
for i = 2:length(dataB_full)
    if maskB_full(i) > maskB_full(i-1)
        startIndex = i;
    elseif maskB_full(i) < maskB_full(i-1)
        sni = max(1,startIndex-5);
        eni = min(length(dataB_full), i+5);
        neighborMean = mean([dataB_full(sni:sni+10);dataB_full(eni:eni+10)]);
        if mod(signalCount,2) == 1
            dataB_full(startIndex:i) = neighborMean*0.5*randn(i-startIndex+1,1)+meanB;
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

subLen = ceil(subLen_full/downsampleRatio);
maxWarp = ceil(maxWarp_full/downsampleRatio);


randRatio = 1/10;
%%%%%%%%%%%%%%%
%%%   DTW   %%%
%%%%%%%%%%%%%%%
tic;
[indicator] = chimp(dataA, dataB, subLen, maxWarp, 20)/sqrt(2*subLen);
toc;
% figure; plot(distancesAB);
N = sum((classificationWindowA(2:end) - classificationWindowA(1:end-1)) == 1);
top20Indices = selectTopN(N,-indicator,subLen);
binaryTop20 = zeros(size(indicator));
binaryTop20(top20Indices) = 1;
impulseTop20 = binaryTop20.*indicator;
indicatorShape = normpdf(1:subLen+1,subLen/2+1, subLen/8);
indicatorShape = indicatorShape/max(indicatorShape);
for i=top20Indices
    startIndex = i-floor(subLen/2);
    startTrim = max(0,0-startIndex + 1);
    endIndex = startIndex + length(indicatorShape)-1;
   binaryTop20(startIndex + startTrim:endIndex) = indicatorShape(startTrim+1:end);
end
% figure; plot(max(0, indicator)/max(indicator)); hold on; plot(classificationWindowA); plot(binaryTop20);

mserror = 1-2*mse(binaryTop20/max(binaryTop20),indicator/max(indicator));
dataANorm = dataA-nanmean(dataA);
dataANorm = dataANorm/nanstd(dataANorm);
dataBNorm = dataB-nanmean(dataB);
dataBNorm = dataBNorm/nanstd(dataBNorm);
complexity = (nanmean(abs(diff(dataANorm))) + nanmean(abs(diff(dataBNorm))))/4;
confidence = mserror*complexity;
% mse(zscore(impulseTop20),zscore(indicator))

%%%%%%%%%%%%%%
%%%   ED   %%%
%%%%%%%%%%%%%%
% numRandTries = length(dataA)*randRatio;
% N = 20;
% distancesAA_ED = nan(length(dataA)-subLen,1);
% for indexA1 = 1:length(distancesAA)
%     if mod(indexA1,100) == 0
%         sprintf("AA: index %d", indexA1)
%     end
%     tempDistances = [];
%     ssA = zscore(dataA(indexA1:indexA1+subLen-1));
%     
%     while length(tempDistances) < numRandTries
%         indexA2 = ceil(rand()*(length(dataA)-subLen));
% %         sprintf("AA: index2 %d", indexA2)
% %         sprintf("index diff: %d", abs(indexA2 - indexA1))
%         if abs(indexA2 - indexA1) > subLen/2
% %            sprintf("AA: index2 %d", indexA2)
%            ssA2 = zscore(dataA(indexA2:indexA2+subLen-1));
%            tempDistances = [tempDistances,norm(ssA-ssA2)];%dtw_upd(ssA, ssA2,0);
%         end
%     end
%     tempDistancesSort = sort(tempDistances);
%     knn = tempDistancesSort(1:N);
%     distancesAA_ED(indexA1) = nanmean(knn);
% end
% % figure; plot(distancesAA);
% 
% 
% 
% numRandTries = length(dataB)*randRatio;
% distancesAB_ED = nan(length(dataA)-subLen,1);
% for indexA1 = 1:length(distancesAB)
%     if mod(indexA1,100) == 0
%         sprintf("AB: index %d", indexA1)
%     end
%     tempDistances = [];
%     ssA = zscore(dataA(indexA1:indexA1+subLen-1));
%     
%     while length(tempDistances) < numRandTries
%         indexB = ceil(rand()*(length(dataB)-subLen));
% %         if abs(indexA2 - indexA1) > subLen/2
%            ssB = zscore(dataB(indexB:indexB+subLen-1));
%            tempDistances = [tempDistances,norm(ssA - ssB)];%dtw_upd(ssA, ssA2,0);
% %         end
%     end
%     tempDistancesSort = sort(tempDistances);
%     knn = tempDistancesSort(1:N);
%     distancesAB_ED(indexA1) = nanmean(knn);
% end
% % figure; plot(distancesAB);
% 
% figure; plot(max(0, distancesAB_ED-distancesAA_ED)); hold on; plot(classificationWindowA*3-1);


% topNDTW = selectTopN(N,-1*(distancesAB - distancesAA),ceil(subLen));
% binaryTopNDTW = zeros(size(maskA));
% binaryTopNDTW(topNDTW) = 1;
% 
% 
% topNED = selectTopN(N,-1*(distancesAB_ED - distancesAA_ED),ceil(subLen));
% binaryTopNED = zeros(size(maskA));
% binaryTopNED(topNED) = 1;

% figure; 
% plot(binaryTopNDTW);
% hold on; 
% plot(binaryTopNED); 
% plot(classificationWindowA); 
% plot(dataA/max(abs(dataA)));