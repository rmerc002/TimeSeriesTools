
subLength = 100;
exclusionLength = subLength;

% ts = tsMale(1:100000);
% ts = randn(10000,1);
% sinIndices = 1:1000:10000-subLength;
% for i = 1:length(sinIndices)
%     si = sinIndices(i);
%    ts(si:si+subLength-1) = sin(linspace(0,2*pi,subLength)); 
% end
tsNonNan = ts;
tsNonNan(isnan(ts)) = nanmean(ts);

tic;
[matrixProfile, matrixProfileIdx, isvalidwindow, motifsIdx, discordsIdx] = mpx_v3(tsNonNan, subLength, subLength, false);
endTime = toc;
fprintf("Matrix Profile completed in %.2f seconds\n",endTime);

tic;
[motifCandidateIndices, d] = allLowestDistanceIndices(matrixProfile, subLength, exclusionLength);
plotSortedBinnedSamples(motifCandidateIndices, ts, [], subLength, 10, 10, "");

% minNum = 1;
% maxNum = min(1000, length(motifCandidateIndices));
% subLenSeries = getSubLenSeries(minNum, maxNum, 30);
frequencySeries = [2,4,8,16,32];
frequencySeries = frequencySeries - 1;
motifCandidateScores = zeros(length(motifCandidateIndices), length(frequencySeries));

for i=1:length(motifCandidateIndices)
    motifCandidateIndex = motifCandidateIndices(i);
    motif = tsNonNan(motifCandidateIndex:motifCandidateIndex + subLength - 1);
    motif(isnan(motif)) = nanmean(motif);
    distProfile = real(MASS_V2(tsNonNan, motif)); %fix to support nan
    
    [motifNNCandidateIndices, nnDist] = allLowestDistanceIndices(distProfile, subLength, exclusionLength);
    
    for numNNIndex = 1:length(frequencySeries)
        numNN = min(frequencySeries(numNNIndex), length(nnDist));
        motifCandidateScores(i, numNNIndex) = mean(real(nnDist(2:numNN+1)))/sqrt(2*subLength); %ignore motif in distance profile
    end
    
end
endTime = toc;
fprintf("MP Frequency completed in %.2f seconds\n",endTime);

% plotDatasetFeatures(motifCandidateScores(1:100,:)')

%%% Now output the best candidate for each number of nearest neighbors
[bestCandidateDistances, minIndices] = min(motifCandidateScores, [],1);
bestCandidateIndicesForEachNumNN = motifCandidateIndices(minIndices);
% 
% 
plotSortedBinnedSamples(bestCandidateIndicesForEachNumNN, ts, [], subLength, 10,10,"");
% 
% motifIndex = 10745;
% motif = tsNonNan(motifIndex:motifIndex+subLength-1);
% distProfile = real(MASS_V2(tsNonNan, motif)); %fix to support nan
% [distProfileCandidateIndices, d] = allLowestDistanceIndices(distProfile, subLength, exclusionLength);
% plotSortedBinnedSamples(distProfileCandidateIndices, tsNonNan, [], subLength, 10, 10, "");

%%% Pan Motif Heatmap Preparation
panMotifFrequencies = nan(length(ts), length(frequencySeries));
for indexFreq = 1:length(frequencySeries)
    
   for indexCandi = 1:length(motifCandidateIndices)
       candidateIndex = motifCandidateIndices(indexCandi);
       startIndex = candidateIndex;
       endIndex = startIndex + subLength - 1;
       panMotifFrequencies(startIndex:endIndex, indexFreq) = motifCandidateScores(indexCandi, indexFreq);
   end
end

dummyIndices = ones(size(panMotifFrequencies));
visualizeMMPAB(ts, panMotifFrequencies, dummyIndices, [], [], [], frequencySeries, frequencySeries, "KNN",""); 


%%% Scatter heatmap
%%% X-axis is time
%%% Y-axis is distance
%%% Color is frequency
colors = [1, 0, 0; 252/255, 136/255, 3/255; 235/255, 192/255, 5/255; 0/255, 174/255, 194/255; 0/255  55/255, 194/255];

figure;
hold on;
for indexFreq = 1:length(frequencySeries)
    color = colors(indexFreq,:);

    x = motifCandidateIndices;
    y =  motifCandidateScores(:, indexFreq);
    scatter(x, y, 20,  color, 'filled');

end
hold off;

