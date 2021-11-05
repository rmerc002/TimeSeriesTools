function [indicator] = chimp(dataA, dataB, subLen, maxWarp, N)
randRatio = 1/10;
%%%%%%%%%%%%%%%
%%%   DTW   %%%
%%%%%%%%%%%%%%%
tic;
numRandTries = length(dataA)*randRatio;

K = randRatio*N*20;%numRandTries*N/(length(dataA)/20);
distancesAA = nan(length(dataA)-subLen,1);
for indexA1 = 1:length(distancesAA)
    if mod(indexA1,1000) == 0
        sprintf("AA: index %d", indexA1)
    end
    tempDistances = [];
    ssA = zscore(dataA(indexA1:indexA1+subLen-1));
    
    while length(tempDistances) < numRandTries
        indexA2 = ceil(rand()*(length(dataA)-subLen));
%         sprintf("AA: index2 %d", indexA2)
%         sprintf("index diff: %d", abs(indexA2 - indexA1))
        if abs(indexA2 - indexA1) > subLen/2
%            sprintf("AA: index2 %d", indexA2)
           ssA2 = zscore(dataA(indexA2:indexA2+subLen-1));
           tempDistances = [tempDistances,dtw_upd(ssA,ssA2,maxWarp)];%dtw_upd(ssA, ssA2,0);
        end
    end
    tempDistancesSort = sort(tempDistances);
    knn = tempDistancesSort(1:N); %This N needs to be related to number of samples 
    distancesAA(indexA1) = nanmean(knn);
end
% figure; plot(distancesAA);



numRandTries = length(dataB)*randRatio;
distancesAB = nan(length(dataA)-subLen,1);
for indexA1 = 1:length(distancesAB)
    if mod(indexA1,1000) == 0
        sprintf("AB: index %d", indexA1)
    end
    tempDistances = [];
    ssA = zscore(dataA(indexA1:indexA1+subLen-1));
    
    while length(tempDistances) < numRandTries
        indexB = ceil(rand()*(length(dataB)-subLen));
%         if abs(indexA2 - indexA1) > subLen/2
           ssB = zscore(dataB(indexB:indexB+subLen-1));
           tempDistances = [tempDistances,dtw_upd(ssA,ssB,maxWarp)];%dtw_upd(ssA, ssA2,0);
%            tempDistances = [tempDistances,norm(ssA-ssB)];%dtw_upd(ssA, ssA2,0);

%         end
    end
    tempDistancesSort = sort(tempDistances);
    knn = tempDistancesSort(1:N);
    distancesAB(indexA1) = nanmean(knn);
end
toc;
% figure; plot(distancesAB);

indicator = max(0, distancesAB-distancesAA);
end

