%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Relative Frequency Matrix Profile   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = data'; %(samples, dimensions)
plotDatasetFeatures(data); 
figure; plot(data);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Relative Frequency Matrix Profile   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classIndex = 4;
m = 200;
maxFreq = 5; %%% Number of nearest neighbors per time index
ts = tsMove01Norm{classIndex};
[mps] = RelativeFrequencyMatrixProfile(ts, ts, m, maxFreq);

KMotifs = 10; %%%Number of motifs to return
motifIndices = NearestNeighborSelection(mps(end,:), m, KMotifs);

for motifIndex = 1:KMotifs
    plotQueryIndexNN(ts, motifIndices(motifIndex), m, maxFreq+1);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
startLen = 3;
endLen = 1000; 
numStep = 5;
plotFeatureSweep(trainClassTS{classIndex}, 3,1000,5);

for classIndex = 1:length(classes)
    plotDataTransforms(classTS{classIndex}(1:10000));
end

%%%Plot Plato with nearestest neighbors in context
numNN = 5;
plotQueryNN(trainClassTS{classIndex}, classPlatos{classIndex}(2,:), numNN);