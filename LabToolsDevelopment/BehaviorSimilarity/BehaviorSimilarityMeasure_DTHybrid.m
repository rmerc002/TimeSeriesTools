function [scores] = BehaviorSimilarityMeasure_DTHybrid(dataA, dataB, dataC, subLength, distanceThreshold)

dataA = reshape(dataA, length(dataA), 1);
dataB = reshape(dataB, length(dataB), 1);
dataC = reshape(dataC, length(dataC), 1);

distanceThreshold = 0.75;

endRegimeAIndex = length(dataA)-subLength+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mpAB, mpiAB] = mpx_v3([dataA; dataB], floor(subLength/2), subLength, false);
mpAB = mpAB(1:endRegimeAIndex);
mpAB = mpAB/sqrt(2*subLength);
mpiAB = mpiAB(1:endRegimeAIndex);
% mpiAB_Count = mpiAB <= length(endRegimeAIndex);

mpActive = mpAB <= distanceThreshold;
mpAB = mpAB(mpActive);
mpiAB = mpiAB(mpActive);

mpiAB_Count = mpiAB <= endRegimeAIndex;

mpLinspace = linspace(1,0,length(mpiAB_Count));
[~,mpAB_sortedOrder] = sort(mpAB);
mpAB_orderedLinspace = zeros(length(mpLinspace),1);
mpAB_orderedLinspace(mpAB_sortedOrder) = mpLinspace;
mpAB_weightedElements = mpAB_orderedLinspace.*mpiAB_Count;
mpAB_weightedMean = mean(mpAB_weightedElements);


scores.AB = max(0,(mpAB_weightedMean-0.25))*4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mpAC, mpiAC] = mpx_v3([dataA; dataC], floor(subLength/2), subLength, false);
mpAC = mpAC(1:endRegimeAIndex);
mpAC = mpAC/sqrt(2*subLength);
mpiAC = mpiAC(1:endRegimeAIndex);
% mpiAC_Count = mpiAC <= length(endRegimeAIndex);

mpActive = mpAC <= distanceThreshold;
mpAC = mpAC(mpActive);
mpiAC = mpiAC(mpActive);

mpiAC_Count = mpiAC <= endRegimeAIndex;

mpLinspace = linspace(1,0,length(mpiAC_Count));
[~,mpAC_sortedOrder] = sort(mpAC);
mpAC_orderedLinspace = zeros(length(mpLinspace),1);
mpAC_orderedLinspace(mpAC_sortedOrder) = mpLinspace;
mpAC_weightedElements = mpAC_orderedLinspace.*mpiAC_Count;
mpAC_weightedMean = mean(mpAC_weightedElements);


scores.AC = max(0,(mpAC_weightedMean-0.25))*4;