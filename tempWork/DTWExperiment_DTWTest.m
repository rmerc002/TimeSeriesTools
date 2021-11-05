



% general testing variables
testIndex = 1;
offset = testIndex;

%%%%%%%%%%%%%%%%%%%%%
%General | pre-process weak labels
processedWeakLabels = weakLabelData - nanmean(weakLabelData);
processedWeakLabels = processedWeakLabels/nanstd(processedWeakLabels);
binaryLabels = processedWeakLabels(testIndex:end)' > 2;

maskOutIndices = [5189, 7463, 9679, 11895, 14145, 16860, 19195, 21389, 23687, 25996, 29926, 32267, 34586, 36895, 39231, 41496, 43824, 46377, 48823, 51165];
for index = 1:length(maskOutIndices)
    binaryLabels(maskOutIndices(index:index+600)) = 0;
end

sampleRate = 30; %30 Hz
dilateTime = 5; %5 seconds before and after
dilateRadius = dilateTime*sampleRate;
wle = binaryLabels;
if sum(binaryDistances(:,1:dilateRadius)) >= 1
   wle(:,1:dilateRadius) = 1; 
end
for index =dilateRadius+1:length(binaryLabels)-dilateRadius - 1
    if sum(binaryLabels(:,index - dilateRadius+1:index + dilateRadius-1)) >= 1
        wle(index) = 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DTW Specific Test %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% DTW | query test set to get distances
distances{1} = inf(1,length(motionData)-subLen+1 - (testIndex));
for kindex = 1:size(motionData,1)-subLen+1 - (testIndex)
    distances{1}(kindex) = dtw_upd(motif1,zscore(motionData(offset + kindex:offset + kindex+subLen-1,1)),maxWarp);
%     if mod(kindex,1000) == 0
%        fprintf("kindex: %d\n",kindex); 
%     end
end
% dtw_distances = distances; %save for plotting
%%%%%%%%%%%%%%%%%%%%%%
%DTW | threshold 
threshold = DTW_motif_distance*2;
binaryDistances{1} = distances{1} < threshold;

%%%%%%%%%%%%%%%%%%%%%%%%
%DTW | Choose Top 11
topNIndices = selectTopN(20, distances, ceil(subLen/2));
binaryTopN{1} = zeros(1,length(distances));
binaryTopN{1}(topNIndices) = 1;


%%%%%%%%%%%%%%%%%%%%%%%
%DTW | plot
figure;
plot(binaryDistances{1} + 1.05)
% plot(binaryTop11 + 1.05);
hold on;
plot(binaryTopN{1} - 1.05);
plot(wle);
hold off;


%%%%%%%%%%%%%%%%%%%%%%%%
%%% ED Specific Test %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% ED | query test set to get distances
[distances_ed,mpi] = mpx_AB(motionData(testIndex:end), motif1,subLen);
distances{2} = distances_ed*sqrt(2*subLen);

% ed_distances = distances; %save for plotting
%%%%%%%%%%%%%%%%%%%%%%
%ED | threshold 
threshold = ED_motif_distance*2;
binaryDistances = distances{2} < threshold;

%%%%%%%%%%%%%%%%%%%%%%%%
%ED | Choose Top 11
topNIndices = selectTopN(20, distances{2}, ceil(subLen/2));
binaryTopN{2} = zeros(1,length(distances));
binaryTopN{2}(topNIndices) = 1;

%%%%%%%%%%%%%%%%%%%%%%%
%ED | plot
figure;
plot(binaryDistances{2} + 1.05)
% plot(binaryTop11 + 1.05);
hold on;
plot(binaryTopN{2} - 1.05);
plot(wle);
hold off;
