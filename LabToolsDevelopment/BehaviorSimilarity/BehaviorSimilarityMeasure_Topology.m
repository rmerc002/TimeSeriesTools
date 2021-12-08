function [scores] = BehaviorSimilarityMeasure_Topology(dataA, dataB, dataC, subLength, distanceThreshold)

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

mpAB_mean = mean(mpiAB_Count);

scores.AB = max(0, (mpAB_mean-0.5))*2;

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

mpAC_mean = mean(mpiAC_Count);

scores.AC = max(0,(mpAC_mean - 0.5))*2;