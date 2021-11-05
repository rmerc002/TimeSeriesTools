function [shapeletCandidateIndices] = getFilteredShapeletCandidates(positiveTS, negativeTS, subLength, groundTruthIndices, figureSavePath)

    showPlot = true;
    if isempty(figureSavePath) || strcmp(figureSavePath,"")
        showPlot = false;
    end

    [MP_AA, MP_AA_Indices] = mpx_shapeletDiscovery(positiveTS,subLength,subLength);
    MP_AA = real(MP_AA);
    [MP_AB, MP_AB_Indices] = mpx_AB_shapeletDiscovery(positiveTS, negativeTS, subLength);
    MP_AB = real(MP_AB);

    minLength = min([length(MP_AA),length(MP_AB)]); %TODO, this is fishy, I think lengths should be the same
    MP_Diff = MP_AB(1:minLength) - MP_AA(1:minLength);

    [shapeletCandidateIndices] = topSubsequenceIndices_V03(MP_Diff, subLength);
    
    shapeletCandidateIndices = shapeletCandidateIndices(1:2);
    
    shapeletCandidateIndices(2) = MP_AA_Indices(shapeletCandidateIndices(1));
    
    
    numCandidates = length(shapeletCandidateIndices);
    
    fileName = "internVars_candidatesVsGroundTruth";
    filePath = fullfile(figureSavePath, fileName + ".mat");
    save(filePath, 'shapeletCandidateIndices', 'groundTruthIndices');

    if showPlot == true
        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41]; 
        blueColor = [0,0.29,0.73];
        grayColor = [0.65,0.65,0.65];
        lightBlueColor = [0.01, 0.83,0.99];

        tsLength = length(positiveTS);
        maxTSLength = max(length(positiveTS),length(negativeTS));
        %%% for debugging
        fig = figure; 
        set(gcf, 'Position', [0,100,2050,400]);
        
        tiledlayout(5,1);
        
        ax1 = nexttile;
%         subplot(5,1,1);
        plot(negativeTS);
        formattedTitle = sprintf("Negative Time Series, Length=%d",length(negativeTS));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        xlim([1,maxTSLength]);

        ax2 = nexttile;
%         subplot(5,1,2);
        plot(positiveTS);
        hold on;
        for i = 1:length(groundTruthIndices)
            gti = groundTruthIndices(i);
            %%%TODO: Change the size back to subLength instead of
            %%%2*subLength
            %%% when ground truth is not concatenated
            plot(gti:gti+subLength-1,positiveTS(gti:gti+subLength-1),'Color',greenColor);
        end
        hold off;
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("Positive Time Series, Length=%d",length(positiveTS));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        xlim([1,maxTSLength]);
        
        ax3 = nexttile;
%         subplot(5,1,3);
        plot(MP_AA,'Color',blueColor);
        hold on;
        plot(MP_AB,'Color',redColor);
        xlim([1,maxTSLength]);
        formattedTitle = "\color{redColor}MP_AB, \color{blueColor}MP_AA";
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);

        ax4 = nexttile;
%         subplot(5,1,4);
        plot(MP_Diff,'Color',grayColor);
        hold on;
        scatter(shapeletCandidateIndices, max(MP_Diff)*ones(1,length(shapeletCandidateIndices)),10,'MarkerFaceColor',lightBlueColor,'MarkerEdgeColor',lightBlueColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
        hold off;
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("MP_Diff");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        ax5 = nexttile;
%         subplot(5,1,5);
        plot(0,0);
        hold on;
        for gti = groundTruthIndices
            plot([gti,gti],[0,1],'Color',greenColor);
        end
        xlim([1,tsLength]);
        formattedTitle = sprintf("GroundTruth");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
        
        
        figureFileName = "getFilteredShapeletCandidates";
        saveas(fig, fullfile(figureSavePath, figureFileName + ".fig"));
        saveas(fig, fullfile(figureSavePath, figureFileName + ".png"));
        close(fig);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Plot Effectiveness of Candidate Selection %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        numCandidatesToShow = min(10, length(shapeletCandidateIndices));
        x = zeros(numCandidatesToShow,2*subLength);
        y = zeros(size(x));

        cmappingGlobal = false(1,length(positiveTS));
        for i = 1:length(groundTruthIndices)
            gti = groundTruthIndices(i);
            %%%TODO: Change the size back to subLength instead of
            %%%2*subLength
            %%% when ground truth is not concatenated
            cmappingGlobal(gti:gti+subLength-1) = true;
        end
        
        fig = figure;
        plot(0,0);
        hold on;
        x = 1:2*subLength;
        percentGroundTruthInCandidate = zeros(numCandidatesToShow,1);
        for i = 1:numCandidatesToShow
            startIndex = shapeletCandidateIndices(i)-ceil(subLength/2);
            leadingWhiteSpace = max(0,1+ -startIndex);% if startIndex == 0, whitespace = 1
            startIndex = max(1, startIndex);
            endIndex = startIndex + subLength - 1 - leadingWhiteSpace;
            endIndex = min(endIndex, length(positiveTS));
            displayLength = endIndex - startIndex + 1;
            
            y = nan(1,2*subLength);
            y(1+leadingWhiteSpace:1+leadingWhiteSpace + displayLength - 1) = positiveTS(startIndex:endIndex);
            cmapping = zeros(1,2*subLength);
            cmapping(1+leadingWhiteSpace:1+leadingWhiteSpace + displayLength - 1) = cmappingGlobal(startIndex:endIndex);
            
            %%% For accuracy reporting, not plotting
            cmappingTrueLength = cmappingGlobal(shapeletCandidateIndices(i): shapeletCandidateIndices(i) + subLength-1);
            percentGroundTruthInCandidate(i) = sum(cmappingTrueLength)/subLength;
            
            tempMin = min(y);
            tempMax = max(y);
            tempRange = max(1e-5, (tempMax-tempMin));

            xfalse = x;
            xtrue = x;

            xfalse(cmapping == true) = nan;
            xtrue(~cmapping == true) = nan;

            yfalse = y;
            ytrue = y;

            yfalse(cmapping == true) = nan;
            ytrue(~cmapping == true) = nan;

            plot(xfalse,1-i+0.95*(yfalse-tempMin)/tempRange,'Color',[0.5,0.5,0.5]);
            plot(xtrue,1-i+0.95*(ytrue-tempMin)/tempRange,'Color',[0,0.8,0.2]);
            gti = find(cmapping,1);
            if gti
                scatter(gti,1-i+0.95*(ytrue(gti)-tempMin)/tempRange,10,'MarkerFaceColor',[0,0.8,0.2], 'MarkerEdgeAlpha',0);
            end
           
        end
        xlim([1,2*subLength]);
        r = rectangle();
        r.Position = [1, 1-i, ceil(size(x,2)/4), numCandidatesToShow];
        r.FaceColor = [1,1,1, 0.45];
        r.LineStyle = '--';
        r.EdgeColor = [0,0,0, 0.45];

        r = rectangle();
        r.Position = [ceil(size(x,2)/4)*3, 1-i, ceil(size(x,2)/4), numCandidatesToShow];
        r.FaceColor = [1,1,1, 0.45];
        r.LineStyle = '--';
        r.EdgeColor = [0,0,0, 0.45];
        hold off;
        
        legendLabels = {};
        legendMarkers = [];
        for i = 1:length(percentGroundTruthInCandidate)
            value = 100* percentGroundTruthInCandidate(i);
            legendLabels{end+1} = sprintf('%.1f%%',value);
            legendMarkers(end+1) = line(nan, nan, 'Linestyle', 'none', 'Marker', '.', 'Color', 'none');
        end
        
        legend(legendMarkers,legendLabels);
        formattedTitle = sprintf("MP-Diff Candidate Selection Performance");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        figureFileName = "candidateSelectionPerformance";
        saveas(fig, fullfile(figureSavePath, figureFileName + ".fig"));
        saveas(fig, fullfile(figureSavePath, figureFileName + ".png"));
        close(fig);

    end
    
     %%% assure var path exists
     varsDirPath = fullfile(figureSavePath, 'Vars');
     if ~isfolder(varsDirPath)
         mkdir(varsDirPath);
     end
    
     %%% internal vars
     dataSetFileName = "getFilteredShapeletCandidates_internalVars";
     internalVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
     save(internalVarsFilePath, 'MP_AA','MP_AB','MP_AA_Indices','MP_AB_Indices','MP_Diff','shapeletCandidateIndices','numCandidates','percentGroundTruthInCandidate');
     
     %%% output vars
     dataSetFileName = "getFilteredShapeletCandidates_output";
     outputVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
     save(outputVarsFilePath, 'shapeletCandidateIndices');
end