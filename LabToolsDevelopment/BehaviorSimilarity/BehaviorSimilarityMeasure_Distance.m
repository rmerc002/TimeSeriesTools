function [scores] = BehaviorSimilarityMeasure_Distance(dataA, dataB, dataC, subLength)

mpAA = mpx_v3(dataA, floor(subLength/2), subLength, false)/sqrt(2*subLength);
mpAA = min(1, mpAA);
scores.AAmean = mean(mpAA);
scores.AAstd = std(mpAA);

mpBB = mpx_v3(dataB, floor(subLength/2), subLength, false)/sqrt(2*subLength);
mpBB = min(1, mpBB);
scores.BBmean = mean(mpBB);
scores.BBstd = std(mpBB);

mpAC = mpx_ABBA_v2(dataA, dataC, subLength)/sqrt(2*subLength);
mpAC = min(1, mpAC);
scores.ACmean = mean(mpAC);
scores.ACstd = std(mpAC);

mpBC = mpx_ABBA_v2(dataB, dataC, subLength)/sqrt(2*subLength);
mpBC = min(1, mpBC);
scores.BCmean = mean(mpBC);
scores.BCstd = std(mpBC);

mpAB = mpx_ABBA_v2(dataA, dataB, subLength)/sqrt(2*subLength);
mpAB = min(1, mpAB);
scores.ABmean = mean(mpAB);
scores.ABstd = std(mpAB);

mpBA = mpx_ABBA_v2(dataB, dataA, subLength)/sqrt(2*subLength);
mpBA = min(1, mpBA);
scores.BAmean = mean(mpBA);
scores.BAstd = std(mpBA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  boxplot %%%
%%%%%%%%%%%%%%%%
maxLength = max([length(mpAA), length(mpBB), length(mpAC), length(mpBC), length(mpAB), length(mpBA)]);
groupedMatrix = nan(maxLength,6);
% {mpAA, mpBB, mpAB, mpBA, mpAC, mpBC}
groupedMatrix(1:length(mpAA),1) = mpAA(:);
groupedMatrix(1:length(mpAB),2) = mpAB(:);
groupedMatrix(1:length(mpAC),3) = mpAC(:);
groupedMatrix(1:length(mpBB),4) = mpBB(:);
groupedMatrix(1:length(mpBA),5) = mpBA(:);
groupedMatrix(1:length(mpBC),6) = mpBC(:);
boxplot(groupedMatrix,'Labels',{'mpAA', 'mpAB', 'mpAC', 'mpBB', 'mpBA', 'mpBC'});
ylim([-0.1, 1.1]);