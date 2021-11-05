% dataA_store = zeros(floor(size(S1Drill,1)/10)+1,(size(S1Drill,2)));
% dataB_store = zeros(floor(size(S2Drill,1)/10)+1,(size(S2Drill,2)));
% maskA_store = zeros(floor(size(S1Drill,1)/10)+1,(size(S1Drill,2)));
% maskB_store = zeros(floor(size(S2Drill,1)/10)+1,(size(S2Drill,2)));
% classificationWindowA_store = zeros(floor(size(S1Drill,1)/10)+1,(size(S1Drill,2)));
% topNNormal_store = zeros(floor(size(S1Drill,1)/10)+1,(size(S1Drill,2)));
% indicator_store = zeros(floor(size(S1Drill,1)/10)+1,(size(S1Drill,2)));
% 
% error_store = zeros(size(S1Drill,2),1);
% complexity_store = zeros(size(S1Drill,2),1);
% confidence_store = zeros(size(S1Drill,2),1);
% accuracy_store = zeros(size(S1Drill,2),1);
% motionColumn_store = zeros(size(S1Drill,2),1);
% mseResults = zeros(245,1);
for index = 138:245
    motionColumn = index%topConf(index+7);%indices(index)
    randEgyptianWarCardGame;

%     dataA_store(1:length(dataA),motionColumn) = dataA; 
    dataB_store(1:length(dataB),motionColumn) = dataB;
%     maskA_store(1:length(maskA),motionColumn) = maskA;
    maskB_store(1:length(maskB),motionColumn) = maskB;
    classificationWindowA_store(1:length(classificationWindowA),motionColumn) = classificationWindowA;
    indicator_store(1:length(indicator),motionColumn) = indicator;


    topNIndices = selectTopN(20,-indicator_store(:,motionColumn), subLen);
%     binaryTopN = zeros(size(indicator_store,1));
    normalTopN = zeros(size(indicator_store,1),1);
%     binaryTopN(topNIndices) = 1;
    % impulseTopN = binaryTopN.*indicator;
    topNMask = false(size(indicator_store,1),1);
    indicatorShape = normpdf(1:subLen+1,subLen/2+1, subLen/8);
    indicatorShape = indicatorShape/max(indicatorShape);
    for i=topNIndices
        startIndex = i-floor(subLen/2);
        startTrim = max(0,0-startIndex + 1);
        endIndex = startIndex + length(indicatorShape)-1;
       normalTopN(startIndex + startTrim:endIndex) = indicatorShape(startTrim+1:end);
       topNMask(startIndex + startTrim:endIndex) = true;
    end
%     

    topNNormal_store(1:length(normalTopN),motionColumn) = normalTopN;%*max(normalTopN);
%     
    
    
    accuracy_store(motionColumn) = sum(classificationWindowA_store(topNIndices,motionColumn));
   
    mserrorSignal = 1-mse(topNNormal_store(topNMask,motionColumn),indicator_store(topNMask,motionColumn));
    mserrorNoise = 1-mse(zeros(sum(~topNMask),1),indicator_store(~topNMask,motionColumn));

    mserror = mserrorSignal*mserrorNoise;
    dataANorm = dataA_store(:,motionColumn)-nanmean(dataA_store(:,motionColumn));
    dataANorm = dataANorm/nanstd(dataANorm);
    dataBNorm = dataB_store(:,motionColumn)-nanmean(dataB_store(:,motionColumn));
    dataBNorm = dataBNorm/nanstd(dataBNorm);
    complexity = 0.25*(nanmean(abs(diff(dataANorm(maskA_store(:,motionColumn))))) + nanmean(abs(diff(dataBNorm(maskB_store(:,motionColumn))))));
    confidence = mserror*complexity;
    complexity_store(motionColumn) = complexity;
    confidence_store(motionColumn) = confidence;
    error_store(motionColumn) = mserror;
end

% for index = 1:245
%    motionColumn = motionColumn_store(index);
%    classificaitonWindowA = classificationWindowA_store(:,motionColumn);
%    topNIndices = selectTopN(20,-indicator_store(:,motionColumn), subLen/2);
%    accuracy_store(motionColumn) = sum(classificationWindowA(topNIndices));
% end