% counter = 0;
% column = 3;%250;
% for index = 1:size(S1Drill)-1
%     if S1Drill{index,column} ~= S1Drill{index+1,column}
%        fprintf("label %d length %d\n",S1Drill{index,column},counter);
%        counter = 0;
%     else
%         counter = counter +1;
%     end
% end
%now decide on sublength based on desired label

trainDataName = "CMU_MotionCapture"



if trainDataName == "S1Drill"
    motionIndex = 1;%69;
    weakLabelIndex = 2;%135;
    k = 28330; %found to be middle of action 10/10 train test;
    motionDataTrain = S1Drill{:,motionIndex};%64};
elseif trainDataName == "S2Drill"
    motionIndex = 1;%69;
    weakLabelIndex = 2;%135;
    k = 27000;
    motionDataTrain = S2Drill{:,motionIndex};%64};
elseif trainDataName == "S3Drill"
    motionIndex = 1;%69;
    weakLabelIndex = 2;%135;
    k = 35500;
    motionDataTrain = S3Drill{:,motionIndex};%64};
elseif trainDataName == "S4Drill"
    motionIndex = 1;%69;
    weakLabelIndex = 2;%135;
    k = 25000;
    motionDataTrain = S4Drill{:,motionIndex};%64};
elseif trainDataName == "hits"
    motionIndex = 4;
    subLen = 70;
    maxWarp = 14;
    k = ceil(length(hits{:,1}));
    motionDataTrain = hits{:,motionIndex};%64};
    weakLabelIndex = 1;%135;
elseif trainDataName == "CMU_MotionCapture"
    k = ceil(size(train,1));
    motionDataTrain = train;%(:,11) - train(:,1);%64};
    weakLabedlIndex = 1;%135;
    subLen = 400;
    maxWarp = 40;
end

[ED_motif_distance,first_ed_motif,second_ed_motif,DTW_motif_distance,first_dtw_motif,second_dtw_motif] = DTWMotifDiscovery(motionDataTrain(1:k), subLen, maxWarp);

motif1 = zscore(motionDataTrain(first_dtw_motif:first_dtw_motif+subLen-1));
motif2 = zscore(motionDataTrain(second_dtw_motif:second_dtw_motif+subLen-1));


