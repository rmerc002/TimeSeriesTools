% general testing variables
testDataName = "CMU_MotionCapture";
testIndex = k+1;
N = 10;
offset = testIndex;
maxWarp = 60;
% motionDataTrain = S1Drill{:,motionIndex};%64};
% motionDataTrain = hits{3000:15000,1};
% 
% motif1 = zscore(motionDataTrain(first_dtw_motif:first_dtw_motif+subLen-1));
% motif2 = zscore(motionDataTrain(second_dtw_motif:second_dtw_motif+subLen-1));


if testDataName == "S1Drill"
    testData = S1Drill;
    N = 10;
    testIndex = k+1;
    maskOutIndices = [4019, 6359, 8869, 12560, 15120, 17680, 20260, 23030, 25730, 28090, 30390, 32740, 35130, 37770, 40270, 42750, 45140, 47630, 50220, 52760];
    motionData = testData{:,motionIndex};%64};
    weakLabelData = testData{:,weakLabelIndex};%135};
elseif testDataName == "S2Drill"
    testData = S2Drill;
    N = 20;
    testIndex = 1;
    maskOutIndices = [5189, 7463, 9679, 11895, 14145, 16860, 19195, 21389, 23687, 25996, 29926, 32267, 34586, 36895, 39231, 41496, 43824, 46377, 48823, 51165];
    motionData = testData{:,motionIndex};%64};
    weakLabelData = testData{:,weakLabelIndex};%135};
elseif testDataName == "hits"
    testData = hits;
    N = 79;
    testIndex = 1;
    maskOutIndices = [];
    motionData = testData{:,motionIndex};%64};
    weakLabelData = testData{:,weakLabelIndex};%135};
elseif trainDataName == "CMU_MotionCapture"
    testData = train;
    N = 6;
    testIndex = 1;
    maxkOutIndices = [];
    motionData = testData;%(:,motionIndex);%64};
    weakLabelData = testData;%(:,weakLabelIndex);%135};
end





distances = zeros(2,length(motionData)-subLen - testIndex + 1);
binaryDistances = zeros(2,length(motionData)-subLen - testIndex+1);
binaryTopN = zeros(2,length(motionData)-subLen - testIndex + 1);
%%%%%%%%%%%%%%%%%%%%%
%General | pre-process weak labels
processedWeakLabels = weakLabelData;
% processedWeakLabels = [processedWeakLabels(subLen:end);zeros(subLen-1,1)];
processedWeakLabels = processedWeakLabels- nanmean(processedWeakLabels);
processedWeakLabels = processedWeakLabels/nanstd(processedWeakLabels);
binaryLabels = abs(processedWeakLabels)' > 2;


for index = 1:length(maskOutIndices)
    moe = maskOutIndices(index);
    binaryLabels(moe:moe+600) = 0;
end

sampleRate = 30; %30 Hz
dilateTime = 2; %2 seconds before and after
dilateRadius = dilateTime*sampleRate;
wle = binaryLabels;
if sum(binaryLabels(:,1:dilateRadius)) >= 1
   wle(:,1:dilateRadius) = 1; 
end
for index =dilateRadius+1:length(binaryLabels)-dilateRadius - 1
    if sum(binaryLabels(:,index - dilateRadius+1:index + dilateRadius-1)) >= 1
        wle(index) = 1;
    end
end
wle(end-subLen+1:end) = [];
if testIndex > 1
   wle(1:testIndex-1) = []; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DTW Specific Test %%%
%%%%%%%%%%%%%%%%%%%%%%%%%

% DTW | query test set to get distances
% distances(1,:) = inf(1,length(motionData)-subLen+1 - (testIndex));
for kindex = 1:size(motionData,1)-subLen+1 - (testIndex)
    distances(1,kindex) = dtw_upd(motif1,zscore(motionData(testIndex-1 + kindex:testIndex-1 + kindex+subLen-1,1)),maxWarp);
%     if mod(kindex,1000) == 0
%        fprintf("kindex: %d\n",kindex); 
%     end
end
% dtw_distances = distances; %save for plotting
%%%%%%%%%%%%%%%%%%%%%%
%DTW | threshold 
threshold = DTW_motif_distance*3;
binaryDistances(1,:) = distances(1,:) < threshold;

%%%%%%%%%%%%%%%%%%%%%%%%
%DTW | Choose Top 11
topNIndices = selectTopN(N, distances(1,:), ceil(subLen/2));
% binaryTopN(1,:) = zeros(1,length(distances(1,:)));
binaryTopN(1,topNIndices) = 1;


%%%%%%%%%%%%%%%%%%%%%%%
%DTW | plot
figure;
hold on;
% plot(bitand(binaryDistances(1,:),wle) + 2.05);
% plot(binaryDistances(1,:) + 1.05)
plot(motionData(testIndex:end)/max(abs(motionData(testIndex:end))) + 1.05);
plot(binaryTopN(1,:) - 1.05);
plot(wle);
plot(bitand(binaryTopN(1,:),wle) - 2.05);
hold off;
title("DTW TEST");

%%%%%%%%%%%%%%%%%%%%%%%%
%%% ED Specific Test %%%
%%%%%%%%%%%%%%%%%%%%%%%%
motif1 = zscore(motionDataTrain(first_dtw_motif:first_dtw_motif+subLen-1));
motif2 = zscore(motionDataTrain(second_dtw_motif:second_dtw_motif+subLen-1));
ED_motif_distance = norm(motif1'-motif2');

%%%%%%%%%%%%%%%%%%%%%%%%%
% ED | query test set to get distances
[distances_ed,mpi] = mpx_AB(motionData(testIndex:end), motif1,subLen);
distances_ed = distances_ed'*sqrt(2*subLen);
distances(2,:) = distances_ed(1:end-1);

% ed_distances = distances; %save for plotting
%%%%%%%%%%%%%%%%%%%%%%
%ED | threshold 
threshold = ED_motif_distance*3;
binaryDistances(2,:) = distances(2,:) < threshold;

%%%%%%%%%%%%%%%%%%%%%%%%
%ED | Choose Top 11
topNIndices = selectTopN(N, distances(2,:), ceil(subLen/2));
binaryTopN(2,:) = zeros(1,length(distances(2,:)));
binaryTopN(2,topNIndices) = 1;

%%%%%%%%%%%%%%%%%%%%%%%
%ED | plot
figure;
hold on;
% plot(bitand(binaryDistances(2,:),wle) + 2.05);
% plot(binaryDistances(2,:) + 1.05)
plot(motionData(testIndex:end)/max(abs(motionData(testIndex:end,:))) + 1.05);
plot(binaryTopN(2,:) - 1.05);
plot(wle);
plot(bitand(binaryTopN(2,:),wle) - 2.05);
hold off;
title("ED TEST");
