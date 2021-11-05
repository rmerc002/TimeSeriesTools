function [queryShapelet] = getCentralSubsequenceOfTopCluster(positiveTS, shapeletCandidateIndices, subLength, numSignals, figureSavePath)
    showPlot = true;
    if isempty(figureSavePath) || strcmp(figureSavePath,"")
        showPlot = false;
    end
    
    %simple version just picks the top 3 peaks (numSignals)

    queryIndex1 = shapeletCandidateIndices(1);
    queryPart1 = positiveTS(queryIndex1:queryIndex1 + subLength - 1);

    queryIndex2 = shapeletCandidateIndices(2);
    queryPart2 = positiveTS(queryIndex2:queryIndex2 + subLength - 1);

    queryShapelet = (zscore(queryPart1) + zscore(queryPart2))/2;
    
    %%% visualizations 
    subsequences = zeros(sum(~isnan(shapeletCandidateIndices)),subLength);
    for i=1:sum(~isnan(shapeletCandidateIndices))
        subsequences(i,:) = positiveTS(shapeletCandidateIndices(i):shapeletCandidateIndices(i)+subLength - 1);
    end

    if showPlot == true
        fig = figure;
        set(gcf, 'Position', [0,100,200,800]);
        % subplot(10,1,[1,9]);
        plot(0,0); hold on;
        inset = 0.9;
        pi = 0;
        for i = 0:length(shapeletCandidateIndices)
            
            
            if i == 0%queryPosIndex
                lineWidth = 5;
                tempTS = queryShapelet;
                plot([1,length(queryShapelet)],[-0.05,-0.05],'--','Color',[0.8,0.8,0.8]);
            else
                lineWidth = 1;
                ti = shapeletCandidateIndices(i);
                tempTS = positiveTS(ti:ti+subLength-1);
            end
            
            tempMin = min(tempTS);
            tempMax = max(tempTS);
            tempRange = max(1e-5, tempMax-tempMin);

            plot(-pi+inset*(tempTS-tempMin)/tempRange,'LineWidth',lineWidth);

            pi = pi + 1;
        end
        hold off;
        xlim([0,subLength]);
        formattedTitle = sprintf("Candidate Subsequences. Query = mean(top2), Length=%d",length(subLength));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);


        figureFileName = "getCentralSubsequenceOfTopCluster";
        saveas(fig, fullfile(figureSavePath, figureFileName + ".fig"));
        saveas(fig, fullfile(figureSavePath, figureFileName + ".png"));
        close(fig);

        % subplot(10,1,10);
        % plot(queryShapelet);
    end
    
    
    %%% assure var path exists
    varsDirPath = fullfile(figureSavePath, 'Vars');
    if ~isfolder(varsDirPath)
     mkdir(varsDirPath);
    end

    %%% I did not think any internal vars were worthwhile
%     %%% internal vars
%     dataSetFileName = "getCentralSubsequenceOfTopCluster_internal";
%     internalVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
%     save(internalVarsFilePath, 'MP_AA','MP_AB','MP_Diff','shapeletCandidateIndices','numCandidates');

    %%% output vars
    dataSetFileName = "getCentralSubsequenceOfTopCluster_output";
    outputVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
    save(outputVarsFilePath, 'queryShapelet');
end