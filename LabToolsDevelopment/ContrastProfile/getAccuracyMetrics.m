function [tpCount,fpCount,fnCount, tpIndices, fpIndices, fnIndices] = getAccuracyMetrics(positiveTS, queryShapelet, MP_AQ, queryThreshold, subLength, groundTruthIndices, figureSavePath)
    showPlot = true;
    disp(groundTruthIndices);
    tpCount = 0;
    fpCount = 0;
    fnCount = 0;

    tpIndices = [];
    fpIndices = [];
    fnIndices = [];

    MP_AQ_Thresh = MP_AQ .* (MP_AQ <= queryThreshold);

    resultIndices = topSubsequenceIndices_V03(MP_AQ_Thresh, subLength);
    resultIndices = resultIndices(MP_AQ_Thresh(resultIndices) > 0);

    groundTruthFound = zeros(length(groundTruthIndices),1);

    for i1 = 1:length(resultIndices)
        ri = resultIndices(i1);
        isResultMatched = false;
       for i2 = 1:length(groundTruthIndices)
           gti = groundTruthIndices(i2);
           if groundTruthFound(i2) == 0 && abs(ri - gti) <= subLength
                   tpCount = tpCount + 1;
                   tpIndices = [tpIndices,ri];
                   groundTruthFound(i2) = 1;
                   isResultMatched = true;
                   break;
           end
           if groundTruthFound(i2) == 1 && abs(ri - gti) <= subLength
                   fpCount = fpCount + 1;
                   fpIndices = [fpIndices,ri];
                   isResultMatched = true;
                   break;
           end
           %%%if groundTruthFound(i2) == 0 && abs(ri - gti) > subLength
           %%% then I don't care

           %%% if groundTruthFound(i2) == 1 && abs(ri - gti) > subLength
           %%% then I don't care
       end
       if isResultMatched == false
           fpCount = fpCount + 1;
           fpIndices = [fpIndices, ri];
       end
    end

    fnIndices = groundTruthIndices((1-groundTruthFound)== 1);
    fnCount = length(fnIndices);
    
    if showPlot == true
        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41]; 
        blueColor = [0,0.29,0.73];
        grayColor = [0.65,0.65,0.65];

        tsLength = length(MP_AQ) + subLength -1;
        %%% Plots
        fig = figure;
        set(gcf, 'Position', [0,100,2050,400]);
        
        tiledlayout(4,1);
        
        %%% Query
        ax1 = nexttile;
%         subplot(4,1,1);
        plot(queryShapelet);
        formattedTitle = sprintf("Query Shapelet, Length=%d",length(queryShapelet));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        ax2 = nexttile;
%         subplot(4,1,2);
        plot(positiveTS);
        xlim([1,length(positiveTS)]);
        formattedTitle = sprintf("Positive Time Series, Length=%d",length(positiveTS));
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        ax3 = nexttile;
%         subplot(4,1,3);
        plot(MP_AQ);
        hold on;
        plot([1,tsLength],[queryThreshold, queryThreshold],'--','Color',grayColor);
        hold off;
        xlim([1,tsLength]);
        formattedTitle = sprintf("Profile Distance: Query and Test Time Series");
        title(formattedTitle,'Interpreter','none');
        set(gca,'xtick',[],'ytick',[]);

        resultTS = zeros(tsLength, 1);
        for ri = resultIndices
           resultTS(ri) = 1;
        end
        
        ax4 = nexttile;
%         subplot(4,1,4);
        plot(0,0);
        hold on;
        %TODO: There's a bug here somewhere
        %Error using horzcat
        %Dimensions of arrays being concatenated are not consistent.
        rectangle('Position',[1,0.5,tsLength,0.5],'FaceColor',redColor,'EdgeColor',redColor)
        for gtii = 1:length(groundTruthIndices)
           gti = groundTruthIndices(gtii);
           startIndex = max(0, gti - subLength);
           endIndex = min(tsLength, gti + subLength);
           width = max(1,endIndex - startIndex + 1);
           rectangle('Position',[startIndex, 0.5, width,0.5],'FaceColor',greenColor,'EdgeColor',greenColor);
        end

        for tpii = 1:length(tpIndices)
           tpi = tpIndices(tpii);
           plot([tpi,tpi],[0,0.5],'Color',greenColor,'LineWidth',1); 
        end
        for fpii = 1:length(fpIndices)
           fpi = fpIndices(fpii);
           plot([fpi,fpi],[0,0.5],'Color',redColor,'LineWidth',1); 
        end
        fprintf("fnIndices: \n");
        disp(fnIndices);
        for fnii = 1:length(fnIndices)
            fni = fnIndices(fnii);
            fprintf("fni index: %d\n",fni);
           plot([fni,fni],[0,0.5],'--','Color',greenColor,'LineWidth',1);
        end
        hold off;
        xlim([1,tsLength]);
        ylim([0,1.25]);
        formattedTitle = sprintf("Ground Truth(top) vs Results(bottom)");
        title(formattedTitle);
        set(gca,'xtick',[],'ytick',[]);
        
        linkaxes([ax2 ax3 ax4],'x') %ax1 is not on the same scale
        
        figureFileName = "getAccuracyMetrics";
        saveas(fig, fullfile(figureSavePath, figureFileName + ".fig"));
        saveas(fig, fullfile(figureSavePath, figureFileName +".png"));
        close(fig);
        
    end
    
    %%% TODO: remove after test 18, redundant
    dataSetFileName = "getAccuracyMetrics_dataVariables";
    dataSetFilePath = fullfile(figureSavePath, dataSetFileName + ".mat");
    save(dataSetFilePath, 'tpIndices','fpIndices','fnIndices');
    
    
    %%% assure var path exists
    varsDirPath = fullfile(figureSavePath, 'Vars');
    if ~isfolder(varsDirPath)
    	mkdir(varsDirPath);
    end

    %%% I did not think any internal vars were worthwhile
%     %%% internal vars
%     dataSetFileName = "getAccuracyMetrics_internal";
%     internalVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
%     save(internalVarsFilePath, 'MP_AQ','groundTruthIndices');

    %%% output vars
    dataSetFileName = "getAccuracyMetrics_output";
    outputVarsFilePath = fullfile(varsDirPath,dataSetFileName + ".mat");
    save(outputVarsFilePath, 'tpIndices','fpIndices','fnIndices');
end