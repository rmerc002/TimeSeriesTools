function [queryThreshold, thresholdConfidence] = getQueryThresholdAndConfidence(testTS, queryShapelet, numSignals, figureSavePath)
    
    showPlot = true;
    if isempty(figureSavePath) || strcmp(figureSavePath,"")
        showPlot = false;
    end
    %%% Temporarily disable while using ideal threshold
    showPlot=false;
    MP_AQ = [];
    topIndices = [];
    numCandidates = [];
    topAmplitudes = [];
    normalizedTopAmplitudes = [];
    
    %%% I know there will be 3, so take the average distance of the 3rd and 4th
    %%% subsequence
    subLength = length(queryShapelet);
    if size(queryShapelet,1) == 1
        queryShapelet = queryShapelet';
    end
    
%     profileDistances = real(mpx_AB_shapeletDiscovery(testTS, queryShapelet, subLength));
%     topIndices = topSubsequenceIndices_V03(-MP_AQ, subLength);
%     numCandidates = length(topIndices);
%     topIndices = topIndices(1:ceil(numCandidates/2));

    
    %%% dividing by sqrt(2*subLength) normalizes the distances between 0 and 1
    %%% It would be better if there was a higher confidence if low variance on
    %%%   to the left of query and also on the right, even if the distances are
    %%%   close

    %%% I can improve it by first normalizing the topPeaks to be between 0 - 1
    %%%   Then compare the normalized distances
%     topAmplitudes = MP_AQ(topIndices);
% 
%     dist1 = topAmplitudes(numSignals);
%     dist2 = topAmplitudes(numSignals + 1);
% 
%     queryThreshold = (dist1 + dist2)/2;

    [idealThreshold] = getIdealThreshold(queryShapelet, testTS, subLength,numSignals, figureSavePath);
    queryThreshold = idealThreshold;
    
    thresholdConfidence = -1; %Ignore this for now

%     tempMin = min(topAmplitudes);
%     tempMax = max(topAmplitudes);
%     tempRange = max(1e-5, tempMax - tempMin);
%     normalizedTopAmplitudes = (topAmplitudes-tempMin)/tempRange;
% 
%     dist1Norm = normalizedTopAmplitudes(numSignals);
%     dist2Norm = normalizedTopAmplitudes(numSignals + 1);
%     thresholdConfidence = (dist2Norm - dist1Norm);

    if showPlot == true
        fig = figure;
        subplot(3,1,1);
        plot(MP_AQ);
        hold on;
        plot([0,length(MP_AQ)], [queryThreshold,queryThreshold],'--','Color',[0.75,0.75,0.75]);
        hold off;
        formattedTitle = sprintf("MP_AQ  = Distances between testTS and query Q, Length=%d",length(MP_AQ));
        title(formattedTitle,'Interpreter','none');
        set(gca,'xtick',[],'ytick',[]);

        subplot(3,1,2);
        plot(topAmplitudes);
        hold on;
        plot([0,numSignals + 0.5], [queryThreshold,queryThreshold],'--','Color',[0.75,0.75,0.75]);
        hold off;
        formattedTitle = sprintf("MP_AQ, %d Top Amplitudes, sorted",length(normalizedTopAmplitudes));
        title(formattedTitle,'Interpreter','none');
        set(gca,'xtick',[],'ytick',[]);

        subplot(3,1,3);
        plot(normalizedTopAmplitudes);
        hold on;
        yValue = thresholdConfidence + topAmplitudes(numSignals);
        %%% horizontal line
        plot([numSignals,numSignals + 1], [topAmplitudes(numSignals),topAmplitudes(numSignals)],'--','Color',[0.65,0.65,0.65]);
        %%% vertical line
        plot([numSignals + 1,numSignals + 1], [yValue,topAmplitudes(numSignals)],'--','Color',[0.65,0.65,0.65]);
        hold off;
        ylim([-0.1, 1.1]);
        title("normalized Top Amplitudes for selecting confidence. Defined over [0,1]");
        set(gca,'xtick',[],'ytick',[]);
        
        figureFileName = "getQueryThresholdAndConfidence";
        saveas(fig, fullfile(figureSavePath, figureFileName + ".fig"));
        saveas(fig, fullfile(figureSavePath, figureFileName +".png"));
        close(fig);
    end

    %%% assure var path exists
    varsDirPath = fullfile(figureSavePath, 'Vars');
    if ~isfolder(varsDirPath)
    	mkdir(varsDirPath);
    end

    %%% I did not think any internal vars were worthwhile
    %%% internal vars
    dataSetFileName = "getQueryThresholdAndConfidence_internal";
    internalVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
    save(internalVarsFilePath, 'MP_AQ','topIndices','numCandidates','topAmplitudes','normalizedTopAmplitudes');

    %%% output vars
    dataSetFileName = "getQueryThresholdAndConfidence_output";
    outputVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
    save(outputVarsFilePath, 'queryThreshold', 'thresholdConfidence');
end