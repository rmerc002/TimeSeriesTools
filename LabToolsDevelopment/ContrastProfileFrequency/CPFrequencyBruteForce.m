
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
K = 100;
PMP = nan(K, length(ts));
PMPi = nan(K, length(ts));
for indexTime = 1:length(ts) - subLength + 1
    subsequence = ts(indexTime:indexTime + subLength - 1);
   distProfile = MASS_V2(ts, subsequence);
   
   startIndex = max(1, indexTime-subLength);
   endIndex = min(length(ts), indexTime + subLength - 1);
   distProfile(startIndex:endIndex) = inf;
   %%%TODO: get candidates so exclusion zones are used
   [sortedIndices, sortedDist] = allLowestDistanceIndices(distProfile, subLength, exclusionLength);
   PMP(:, indexTime) = sortedDist(1:K);
   PMPi(:, indexTime) = sortedIndices(1:K);
   disp(indexTime);
end
PMPMean = nan(size(PMP));
for indexK = 1:K
   PMPMean(indexK, :) = mean(PMP(1:indexK,:), 1); 
end
endTime = toc;
fprintf("MP Frequency completed in %.2f seconds\n",endTime);


visualizeMMPAB(ts, PMPMean, PMPi, [], [], [], 100, 1:100, "KNN","");
% plotDatasetFeatures(motifCandidateScores(1:100,:)')

%%% Now output the best candidate for each number of nearest neighbors
[bestCandidateDistances, minIndices] = min(PMPMean, [],2);
% 
% 
plotSortedBinnedSamples(minIndices, ts, [], subLength, 10,10,"");
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


