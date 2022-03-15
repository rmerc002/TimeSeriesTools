%%% Input:
%%%   obj.positiveTS: A time series containing at least two instances of a desired behavior
%%%   obj.negativeTS: A time series containing zero instances of a desired behavior
%%%   m: The approximate window size of the desired behavior
%%%   forcePlot: (optional) Display the following two plots
%%%
%%% Output:
%%%   novelet: The subsequence that most distinguishes 
%%%          obj.positiveTS from obj.negativeTS
%%%   obj.noveletIndices: The starting index of each of the K Novelets
%%%   noveletPrimaryContrast: Contrast value of each novelet in the K=1 obj.NP
%%%   noveletNaryContrast: Contrast value of each novelet after appending
%%%     the previous Novelet to obj.negativeTS. May be helpful in identifying
%%%     diminishing returns and retundant behaviors.
classdef OnlineNovelProfile
    properties
        positiveTS
        negativeTS
        mm % Subsequence length
        startLocalOffset
        endLocalOffset

        exclusionLength
        contextLength
        noveltyThreshold
        lMP_AA
        rMP_AA
        lMP_AA_Indices
        rMP_AA_Indices
        MP_AB
        CP % Contrast Profile
        preNP %Highlights seconds instance of behavior
        NP % Novelet Profile
        rNP

        minBufferLength
        oldNewsAnnotation
        firstFoundAnnotation
        secondFoundAnnotation
        bufferStartIndex

        novelets
        noveletIndices
        noveletNNIndices
        noveletScores
    end
    methods
        function obj = OnlineNovelProfile(positiveTS, mm, negativeTS, noveltyThreshold, startLocalOffset, endLocalOffset)
            if nargin < 1
               positiveTS = [];
            end
            if nargin < 2
                mm = 10;
            end
            if nargin < 3 || length(negativeTS) < mm
                negativeTS = [];
            end
            if nargin < 4
                noveltyThreshold = 0;
            end
            if nargin < 5
                startLocalOffset = 0;
            end
            if nargin < 6
                startLocalOffset = 0;
                endLocalOffset = 0;
            end
            
%             obj.positiveTS = reshape(positiveTS, length(positiveTS), 1);
            obj.mm = mm;
            obj.startLocalOffset = startLocalOffset;
            obj.endLocalOffset = endLocalOffset;

            obj.exclusionLength = ceil(obj.mm);
            obj.contextLength = ceil(obj.mm/2);
            obj.negativeTS = reshape(negativeTS, length(negativeTS), 1);
            obj.noveltyThreshold = noveltyThreshold;
            
            obj.minBufferLength = obj.mm*4;
            obj.bufferStartIndex = 1;
            obj.oldNewsAnnotation = false(length(obj.positiveTS),1); %%% boolean
            obj.novelets = []; %%% dimension: (numNovelets, mm)
            obj.noveletIndices = [];
            obj.noveletNNIndices = [];
            obj.noveletScores = [];

            if isempty(obj.negativeTS)
                obj.negativeTS = rand(obj.mm+1,1);
            end

            if noveltyThreshold == 0
                [~,~,cp] = ContrastProfile(positiveTS', obj.negativeTS, obj.mm, true);
                [~, peakValues] = exclusionZonePeaks(cp, obj.mm, obj.exclusionLength);
                figure;
                histogram(peakValues,100);
                userMessage = "";
                userMessage = userMessage + "Please decide on a novelty threshold\n";
                obj.noveltyThreshold = input(userMessage);
            end
            
            obj.bufferStartIndex = 1;

            obj = obj.update(positiveTS);

            end

        function obj = discoverNovelets(obj, t, MP_AB, MP_AA, MP_AA_Indices)
            t = reshape(t, length(t), 1);

            
            bufferLength = length(t);
            initialBufferStartIndex = obj.bufferStartIndex;
            while bufferLength >= obj.minBufferLength
                bufferIndex = obj.bufferStartIndex - initialBufferStartIndex + 1;
                startIndexActionWindow = bufferIndex;
                endIndexActionWindow = min(length(obj.positiveTS), startIndexActionWindow + obj.minBufferLength - 1);
                endIndexMPActionWindow = endIndexActionWindow - obj.mm + 1;
                tSub = t(startIndexActionWindow:endIndexActionWindow);
                MP_ABSubOrig = MP_AB(startIndexActionWindow:endIndexMPActionWindow);
                MP_ABSubUpdated = MP_ABSubOrig;
                MP_AASub = MP_AA(startIndexActionWindow:endIndexMPActionWindow);
                MP_AA_IndicesSub = MP_AA_Indices(startIndexActionWindow:endIndexMPActionWindow);

                %%% Update the MP_AB here
                try
                    for jj = 1:size(obj.novelets,1)
                        novelet = obj.novelets(jj, :)';
                        newMP_AB = mpx_ABBA_v2(tSub, novelet, obj.mm); %%%TODO: verify if this overlap is correct
                        newMP_AB = clipMatrixProfileAmplitude(newMP_AB, obj.mm); 
                        MP_ABSubUpdated(1:length(newMP_AB)) = min(MP_ABSubUpdated(1:length(newMP_AB)), newMP_AB);
                    end
                catch
                    return
                end

                %%%Ignore time series near by
                if obj.bufferStartIndex > obj.startLocalOffset && obj.startLocalOffset > 0
                    startLocalOffset = obj.bufferStartIndex - obj.startLocalOffset;
                    endLocalOffset = obj.bufferStartIndex - obj.endLocalOffset;
                    localTS = obj.positiveTS(startLocalOffset: endLocalOffset);
                    newMP_AB = mpx_ABBA_v2(tSub, localTS, obj.mm); %%%TODO: verify if this overlap is correct
                    newMP_AB = clipMatrixProfileAmplitude(newMP_AB, obj.mm); 
                    MP_ABSubUpdated(1:length(newMP_AB)) = min(MP_ABSubUpdated(1:length(newMP_AB)), newMP_AB);
                end
                
                CPSub = MP_ABSubOrig - MP_AASub;
                CPSub = normalizeContrastProfileAmplitude(CPSub, obj.mm);

                preNPSub = MP_ABSubUpdated - MP_AASub;
                preNPSub = normalizeContrastProfileAmplitude(preNPSub, obj.mm);
                if any(preNPSub>0)
                    %%%This is an important note for the paper
                    %%% There is the problem of selecting a local optimum
                    %%% peak if the greatest value is on the far right of
                    %%% the current window. The issue is that there may be
                    %%% an even greater value within a subsequence length.
                    %%% This is solved by using an exclusion zone of mm
                    %%% along with min buffer size of 4*mm and selecting
                    %%% the top two peaks within the buffer. By choosing
                    %%% the smaller of the two peaks, which is at least mm 
                    %%% due to the exclusion zone, a possibly more optimal 
                    %%% peak is saved for the next round. If the highest
                    %%% peak is on the far right, and there is an
                    %%% undiscovered higher peak within a subsequence
                    %%% length of mm, then in the next selection cycle, the
                    %%% current highest peak will be pushed away due to the
                    %%% exclusion zone of the next better peak.
                   [peakIndices, peakValues] = exclusionZonePeaks(preNPSub, obj.mm, obj.exclusionLength, 2);
                    peakIndex = peakIndices(1);
                else
                    peakIndex = 1;
                end
                
                peakMP_AB = MP_ABSubUpdated(peakIndex);
                peakMP_AA = MP_AASub(peakIndex);
                peakMP_AA_Index = MP_AA_IndicesSub(peakIndex);

                
                peakContrast = peakMP_AB - peakMP_AA;
                peakContrast = normalizeContrastProfileAmplitude(peakContrast, obj.mm);

                if peakContrast > obj.noveltyThreshold
                    obj.noveletNNIndices(end+1) = obj.bufferStartIndex + peakIndex - 1; %%% Swapped because we find the second instance
                    obj.noveletIndices(end+1) = peakMP_AA_Index;                        %%%   then look back to the first
                    novelet = obj.getSubsequenceWithContext(peakMP_AA_Index, obj.contextLength);
                    obj.novelets(end+1,:) = novelet;

%                     newMP_AB = mpx_ABBA_v2(tSub, novelet, obj.mm); %%%TODO: verify if this overlap is correct
%                     newMP_AB = clipMatrixProfileAmplitude(newMP_AB, obj.mm); 
%                     MP_ABSubUpdated(1:length(newMP_AB)) = min(MP_ABSubUpdated(1:length(newMP_AB)), newMP_AB);
                    
                    obj.noveletScores(end+1) = peakContrast;

                    startIndex = obj.bufferStartIndex + peakIndex - 1;
                    endIndex = startIndex + obj.mm - 1;
                    obj.oldNewsAnnotation(startIndex: endIndex) = true;

                    startIndex = peakMP_AA_Index;
                    endIndex = startIndex + obj.mm - 1;
                    obj.oldNewsAnnotation(startIndex: endIndex) = true;
                end

                
                activeEndIndex = obj.bufferStartIndex + length(MP_AASub)-1;
                obj.lMP_AA(obj.bufferStartIndex:activeEndIndex) = MP_AASub;
                obj.lMP_AA_Indices(obj.bufferStartIndex:activeEndIndex) = MP_AA_IndicesSub;
                obj.MP_AB(obj.bufferStartIndex:activeEndIndex) = MP_ABSubUpdated;
                obj.CP(obj.bufferStartIndex:activeEndIndex) = CPSub;
                obj.preNP(obj.bufferStartIndex:activeEndIndex) = preNPSub;

                %%% try to remove impact of the second instance
%                 if peakContrast > obj.noveltyThreshold
%                     startIndexActionWindow = max(1, obj.noveletNNIndices(end) - obj.contextLength);
%                     endIndexActionWindow = min(length(obj.positiveTS), startIndexActionWindow + obj.minBufferLength + obj.contextLength - 1);
%                     tSub = t(startIndexActionWindow:endIndexActionWindow);
%                     newMP_AB = mpx_ABBA_v2(tSub, novelet, obj.mm); %%%TODO: verify if this overlap is correct
%                     newMP_AB = clipMatrixProfileAmplitude(newMP_AB, obj.mm); 
%                     tempAB = obj.MP_AB(startIndexActionWindow:endIndexActionWindow-obj.mm+1);
%                     tempAB = min(tempAB, newMP_AB);
% 
%                     obj.MP_AB(startIndexActionWindow:endIndexActionWindow-obj.mm+1) = tempAB;
%                 end

                
                obj.bufferStartIndex = obj.bufferStartIndex + peakIndex + obj.exclusionLength;
                bufferLength = length(t) - (obj.bufferStartIndex - initialBufferStartIndex + 1) + 1;
            end

            %%% NP caculated from CP and MP_AA_Indices
            for ii = initialBufferStartIndex: length(obj.NP)
                if ~isnan(obj.lMP_AA_Indices(ii))
                    obj.NP(obj.lMP_AA_Indices(ii)) = max(obj.NP(obj.lMP_AA_Indices(ii)), obj.preNP(ii));
                end
            end

            obj.rNP = normalizeContrastProfileAmplitude(obj.MP_AB - obj.rMP_AA, obj.mm);
            
%             for ii = initialBufferStartIndex:activeEndIndex
%                 if obj.CP(ii) > 0.1 && obj.preNP(ii) < obj.noveltyThreshold
%                     obj.oldNewsAnnotation(ii:ii+obj.mm-1) = true;
%                 end
%             end
            %%%When the contrast cannot possibly rise above the threshold
            for ii = initialBufferStartIndex:length(obj.MP_AB)
                if obj.MP_AB(ii) < obj.noveltyThreshold*sqrt(2*obj.mm)
                    obj.oldNewsAnnotation(ii:ii+obj.mm-1) = true;
                end
            end
            

        end

        function obj = update(obj, t)
            if isempty(t)
                return;
            else
                t = reshape(t, length(t), 1);
                obj.positiveTS = [obj.positiveTS; t];
                obj.lMP_AA = [obj.lMP_AA;2*sqrt(obj.mm)*ones(length(t),1)];
                obj.lMP_AA_Indices = [obj.lMP_AA_Indices;nan(length(t),1)];
                obj.rMP_AA = [obj.rMP_AA;2*sqrt(obj.mm)*ones(length(t),1)];
                obj.rMP_AA_Indices = [obj.rMP_AA_Indices;nan(length(t),1)];
                obj.MP_AB = [obj.MP_AB;2*sqrt(obj.mm)*ones(length(t),1)];
                obj.CP = [obj.CP;zeros(length(t),1)];
                obj.NP = [obj.NP;zeros(length(t),1)];
                obj.oldNewsAnnotation = [obj.oldNewsAnnotation; false(length(t),1)];
            end
        
%             numPossibleSubsequences = length(obj.positiveTS) - obj.mm + 1;
            %%% Assume length(obj.MP_AB) includes m-1 points after last valid value
            currentBufferLength = length(obj.positiveTS) - obj.bufferStartIndex + 1;
            if currentBufferLength < obj.minBufferLength
                numSamplesNeededForUpdate = obj.minBufferLength - currentBufferLength;
%                 fprintf("New samples stored, but %d more samples required to act on buffer.\n", numSamplesNeededForUpdate)
                return;
            end
        
            obj.bufferStartIndex = max(1, obj.bufferStartIndex - 2*obj.exclusionLength);%%%TODO: Temp to debug, but maybe necessary
            unprocessedPositiveTS = obj.positiveTS(obj.bufferStartIndex:end);
            unprocessedLength = length(unprocessedPositiveTS);

            numSubsequences = unprocessedLength - obj.mm + 1;

            newMP_AB = sqrt(2*obj.mm)*ones(unprocessedLength,1);
            if length(obj.negativeTS) >= obj.mm
                tempMP_AB = mpx_ABBA_v2(unprocessedPositiveTS, obj.negativeTS, obj.mm);
            else
                tempMP_AB = [];
            end
            newMP_AB(1:length(tempMP_AB)) = tempMP_AB;
            %%% Start Novelet Comparisons
            numNovelets = size(obj.novelets,1);
            noveletDP = sqrt(2*obj.mm)*ones(numNovelets,unprocessedLength);
            
            if numNovelets > 0
                for noveletIndex = 1:numNovelets
                    novelet = reshape(obj.novelets(noveletIndex,:), obj.mm + 2*obj.contextLength, 1);
                    tempNoveletDist = real(MASS_V2(unprocessedPositiveTS, novelet));
                    tempNoveletDist = clipMatrixProfileAmplitude(tempNoveletDist, obj.mm);
                    noveletDP(noveletIndex,1:length(tempNoveletDist)) = tempNoveletDist;
                end
                [noveletMinDist, ~] = min(noveletDP,[],1);
                noveletMinDist = reshape(noveletMinDist, unprocessedLength, 1);
    
                newMP_AB = min(newMP_AB, noveletMinDist);
            end


            %%%Extending the left only matrix profile can be done by
            %%%combining an AB-join between unprocessed posTS with rest of
            %%%posTS , then the left only self-join of the unprocessed
            %%%posTS
            
            newMP_AA = sqrt(2*obj.mm)*ones(unprocessedLength, 1);
            newMP_AA_Indices = ones(unprocessedLength, 1);

%             [~, tempLeftAA, ~, tempLeftAAIndices] = mpxLeftRight(obj.positiveTS, ceil(obj.mm/2), obj.mm);
%             tempLeftAA = tempLeftAA(obj.bufferStartIndex:end);
%             newMP_AA(1:length(tempLeftAA)) = tempLeftAA;
            tempAB = nan(unprocessedLength,1);
            if obj.bufferStartIndex > obj.mm
                [tempAB, ~, tempABIndices] = mpx_ABBA_v2(unprocessedPositiveTS, obj.positiveTS(1:obj.bufferStartIndex-1), obj.mm);
            end
            [~, tempLeftAA, ~, tempLeftAAIndices] = mpxLeftRight(unprocessedPositiveTS, obj.exclusionLength, obj.mm);
            tempLeftAAIndices = tempLeftAAIndices + obj.bufferStartIndex -1;
            
            for ii = 1:numSubsequences
                if isnan(tempAB(ii))
                    newMP_AA(ii) = tempLeftAA(ii);
                    newMP_AA_Indices(ii) = tempLeftAAIndices(ii);
                elseif isnan(tempLeftAA(ii))
                    newMP_AA(ii) = tempAB(ii);
                    newMP_AA_Indices(ii) = tempABIndices(ii);
                elseif tempAB(ii) < tempLeftAA(ii) 
                    newMP_AA(ii) = tempAB(ii);
                    newMP_AA_Indices(ii) = tempABIndices(ii);
                else %%% if AB is nan or 
                    newMP_AA(ii) =  tempLeftAA(ii);
                    newMP_AA_Indices(ii) = tempLeftAAIndices(ii);
                end
            end
%             newMP_AA(1:obj.mm) = 2*sqrt(obj.mm);
%             obj.bufferStartIndex = obj.bufferStartIndex + 2*obj.exclusionLength;
%             obj = discoverNovelets(obj, unprocessedPositiveTS(2*obj.exclusionLength+1:end), newMP_AB(2*obj.exclusionLength+1:end), newMP_AA(2*obj.exclusionLength+1:end), newMP_AA_Indices(2*obj.exclusionLength+1:end));
            
            
 

            %%%Discover novelets
            obj = discoverNovelets(obj, unprocessedPositiveTS, newMP_AB, newMP_AA, newMP_AA_Indices);
            %%%end of update
        end

        function subsequenceWithContext = getSubsequenceWithContext(obj, subsequenceIndex, contextLength)
            startIndex = subsequenceIndex - contextLength; %%% need to handle when exclusion length is at the start or end
            endIndex = subsequenceIndex + obj.mm - 1 + contextLength;

            indicesRaw = startIndex:endIndex;
            indices = indicesRaw - 1;

            indices = mod(indices, length(obj.positiveTS));
            indices = indices + 1;
    
            subsequenceWithContext = obj.positiveTS(indices);
            
            subsequenceWithContext(indicesRaw <= 0) = nan;
            subsequenceWithContext(indicesRaw > length(obj.positiveTS)) = nan;

        end
        function plot(obj)
            %%%%%%%%%%%%%
            %%% PLOTS %%%
            %%%%%%%%%%%%%    
            obj.oldNewsAnnotation = false(length(obj.positiveTS),1); %%%TODO: remove after supporting in main functions

            redColor = [0.73,0.05,0];
%             greenColor = [0,0.73,0.41]; 
            blueColor = [0,0.29,0.73];
            grayColor = [0.75,0.75,0.75];
            lightGrayColor = [0.9, 0.9, 0.9];
%             lightBlueColor = [0.01, 0.83,0.99];
            noveletColor = [143/255, 226/255, 227/255];
            noveletLabelColor = [56/255, 226/255, 235/255];
            noveletNNColor = [227/255, 143/255, 219/255];
            noveletNNLabelColor = [235/255, 56/255, 217/255];
    
%             tsLength = length(obj.positiveTS);
            maxTSLength = max(length(obj.positiveTS),length(obj.negativeTS));
    
            figure('Name','Emergence Profile: Panel','NumberTitle','off'); 
            set(gcf, 'Position', [0,100,2000,600]);
            
            tiledlayout(4,1);
            
            ax1 = nexttile;
            if ~isempty(obj.negativeTS)
                plot(obj.negativeTS,'Color',redColor);
                formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(-)}\\color{black}: Negative Time Series",redColor(1), redColor(2), redColor(3));
                title(formattedTitle);
                set(gca,'xtick',[1,length(obj.negativeTS)],'ytick',[], 'TickDir','out');
            end
            xlim([1,maxTSLength]);
            box off;
            
            ax2 = nexttile;
            hold on;
            %%% Unmatched
            tempPos = obj.positiveTS;
            tempPos(obj.oldNewsAnnotation) = nan;
%             plot(tempPos, 'LineWidth',1);
            plot(tempPos);

    
            %%% Duplicate and discarded
            tempPos = obj.positiveTS;
            tempPos(~obj.oldNewsAnnotation) = nan;
            plot(tempPos,'Color',lightGrayColor,'LineWidth',1);

            %%% Matched Second instance
            for ii = 1:length(obj.noveletIndices)
                startIndex = obj.noveletNNIndices(ii);
                endIndex = startIndex + obj.mm;
                subsequence = obj.positiveTS(startIndex:endIndex);
                plot(startIndex:endIndex, subsequence,'Color',noveletNNColor,'LineWidth',1)
            end
            
            %%% Matched First instance
            for ii = 1:length(obj.noveletIndices)
                startIndex = obj.noveletIndices(ii);
                endIndex = startIndex + obj.mm;
                subsequence = obj.positiveTS(startIndex:endIndex);
                plot(startIndex:endIndex, subsequence,'Color',noveletColor,'LineWidth',1)
            end
            
            markerHeight = 1.2*max(obj.positiveTS);
            scatter(obj.noveletIndices, markerHeight*ones(length(obj.noveletScores),1) ,20,'v','MarkerFaceColor',noveletLabelColor,'MarkerEdgeColor',noveletLabelColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
            scatter(obj.noveletNNIndices, markerHeight*ones(length(obj.noveletScores),1), 20,'v','MarkerFaceColor',noveletNNLabelColor,'MarkerEdgeColor',noveletNNLabelColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 

            
    
            hold off;
            xlim([1,maxTSLength]);
            formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(+)}\\color{black}: Positive Time Series",blueColor(1),blueColor(2),blueColor(3));
            title(formattedTitle);
            set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[], 'TickDir','out');
            xlim([1,maxTSLength]);
            box off;
            
            ax3 = nexttile;
            plot(obj.lMP_AA,'Color',blueColor);
            hold on;
            plot(obj.MP_AB,'Color',redColor);
            xlim([1,maxTSLength]);
            ylim([0,1.1*sqrt(2*obj.mm)]);
            formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}MP^{(+ -)} AB-join, \\color[rgb]{%f,%f,%f}MP^{(+ +)} Left Self-Join",redColor(1), redColor(2), redColor(3),blueColor(1),blueColor(2),blueColor(3));
            title(formattedTitle);
            set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[0,sqrt(2*obj.mm)], 'TickDir','out');
            box off;
            
            ax4 = nexttile;
            hold on;
            plot([0,maxTSLength],[obj.noveltyThreshold, obj.noveltyThreshold],'--','Color',grayColor);
            plot(obj.preNP,'Color',grayColor);
            scatter(obj.noveletIndices, ones(length(obj.noveletScores),1) ,20,'v','MarkerFaceColor',noveletLabelColor,'MarkerEdgeColor',noveletLabelColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
            scatter(obj.noveletNNIndices, ones(length(obj.noveletScores),1), 20,'v','MarkerFaceColor',noveletNNLabelColor,'MarkerEdgeColor',noveletNNLabelColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
            hold off;
            xlim([1,maxTSLength]);
            ylim([-0.1,1.1]);
            formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Novel Profile, \\color[rgb]{%f,%f,%f}First, \\color[rgb]{%f,%f,%f}Second",grayColor(1), grayColor(2), grayColor(3), noveletLabelColor(1), noveletLabelColor(2), noveletLabelColor(3), noveletNNLabelColor(1), noveletNNLabelColor(2), noveletNNLabelColor(3));
            title(formattedTitle);
            set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[0,1], 'TickDir','out');
            box off;
            
            
            linkaxes([ax1 ax2 ax3 ax4],'x')
            
            %%%Generating unique colors for classes when there can be many classes
            %%% I will assume no more than 1000 results
            numColors = 1000;
            colors = lines(numColors);
            
            figure;
            hold on;
            numNovelets = size(obj.novelets,1);
            startIndex = obj.contextLength+1;
            endIndex = startIndex + obj.mm - 1;
            for ki = 1:numNovelets
                startIndexNN = obj.noveletNNIndices(ki);
                endIndexNN = startIndexNN + obj.mm - 1;
                tempTS = obj.positiveTS(startIndexNN: endIndexNN);
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
                plot(1:obj.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color', 1-0.3*(1-colors(ki,:)));


                tempTS = obj.novelets(ki,startIndex:endIndex);
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
    
                tempContextTS = obj.novelets(ki,:);
                plot(-obj.contextLength+1:obj.mm+obj.contextLength,-ki + 0.9*(tempContextTS - tempMin)/tempRange,'Color',lightGrayColor);
                plot(1:obj.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));
            end
            hold off;
            formattedTitle = sprintf("%d Novelets > %.2f", numNovelets, obj.noveltyThreshold);
            title(formattedTitle);
            set(gca,'xtick',[1,obj.mm],'ytick',[], 'TickDir','out');
            xlim([-obj.contextLength+1,obj.mm+obj.contextLength]);
            box off;

            
            %%%novelets sorted by score
%             [~,sortedIndices] = sort(obj.noveletScores,"descend");
%             figure;
%             hold on;
%             numNovelets = size(obj.novelets,1);
%             startIndex = obj.contextLength+1;
%             endIndex = startIndex + obj.mm - 1;
%             for ki = 1:numNovelets
%                 ii = sortedIndices(ki);
%                 startIndexNN = obj.noveletNNIndices(ii);
%                 endIndexNN = startIndexNN + obj.mm - 1;
%                 tempTS = obj.positiveTS(startIndexNN: endIndexNN);
%                 tempMin = min(tempTS);
%                 tempMax = max(tempTS);
%                 tempRange = tempMax-tempMin;
%                 plot(1:obj.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color', 1-0.3*(1-colors(ki,:)));
% 
% 
%                 tempTS = obj.novelets(ii,startIndex:endIndex);
%                 tempMin = min(tempTS);
%                 tempMax = max(tempTS);
%                 tempRange = tempMax-tempMin;
%     
%                 tempContextTS = obj.novelets(ii,:);
%                 plot(-obj.contextLength+1:obj.mm+obj.contextLength,-ki + 0.9*(tempContextTS - tempMin)/tempRange,'Color',lightGrayColor);
%                 plot(1:obj.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));
%             end
%             hold off;
%             formattedTitle = sprintf("%d Emerging Behaviors > %.2f", numNovelets, obj.noveltyThreshold);
%             title(formattedTitle);
%             set(gca,'xtick',[1,obj.mm],'ytick',[], 'TickDir','out');
%             xlim([-obj.contextLength+1,obj.mm+obj.contextLength]);
%             box off;
            %%% end of plot
        end

        function noveletsWithoutContext = getNoveletsWithoutContext(obj)
            noveletsWithoutContext = obj.novelets(:,obj.contextLength+1:obj.contextLength+obj.mm);
        end
    end
end

function [mp] = clipMatrixProfileAmplitude(mp,m)
    mp = min(sqrt(2*m),mp,'includenan'); %clip the anti-correlated values
    mp = real(mp);
end

function [cp] = normalizeContrastProfileAmplitude(cp, m)
    cp = cp/(sqrt(2*m)); %normalize between 0,1
    %discard negative values. These occur when a subsequence in T+ has a
    %closer nearest neighbor in T-. It's not a behavior of interest for
    %this research.
    cp = max(0, cp,'omitnan'); 
end


function [mp_a, mp_b, mpi_a, mpi_b] = mpx_ABBA_v2(a, b, w)
% (mpx_ABBA_v2 titled from Kaveh Kamgar, 2020-12-24 16:56)
% Code and update formulas are by Kaveh Kamgar. 
% GUI and top k motif critera are based on but not identical to some code
% by Michael Yeh. Implementation details are specified above.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product


% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if isrow(a)
    a = transpose(a);
elseif ~iscolumn(a)
    error('Expected a 1D input for input array a');
end

if isrow(b)
    b = transpose(b);
elseif ~iscolumn(b)
    error('Expected a 1D input for input array b');
end

% matrix profile using cross correlation,
% updated to use the function muinvn, which is a nicer implementation of my
% take on Ogita's work

% save original lengths
subcount_a = length(a) - w + 1;
subcount_b = length(b) - w + 1;

% slightly misleading name. This both identifies valid windows and 
% perturbs flat data with a bit of noise, scaled appropriately. If it flips
% any nearest neighbors, they were already unreliable matches.
[a, isvalid_a, first_a, last_a] = find_valid_windows(a, w);
[b, isvalid_b, first_b, last_b] = find_valid_windows(b, w);

mu_a = moving_mean(a, w);
mu_b = moving_mean(b, w);

invn_a = NaN(length(a) - w + 1, 1);
invn_b = NaN(length(b) - w + 1, 1);

% These are truncated so that they start and end on valid subsequences.
% Without that, any adjustment for missing data is susceptible to invalid
% or perturbed initial comparisons
for i = 1 : length(invn_a)
    if isvalid_a(i)
        invn_a(i) = 1 / norm(a(i : i + w - 1) - mu_a(i)); 
    end
end

for i = 1 : length(invn_b)
    if isvalid_b(i)
        invn_b(i) = 1 / norm(b(i : i + w - 1) - mu_b(i));
    end
end

% This means that for some subsequence, 1/(norm(subseq - mean(subseq)) 
% is not a finite value, AND this subsequence was not noticed by a mising
% data or constant window check. It is intended as an indicator that you
% should investigate the data itself, because we can't concretely identify
% the issue here. It can be something as simple as intermediate overflow,
% but with double precision inputs, even that is typically indicative 
% of a problem with the input. 

if any(~isfinite(invn_a) & isvalid_a)
    warning('Times series a has non-normalizable subsequences, which were not detected by the missing data passes');
end

if any(~isfinite(invn_b) & isvalid_b)
    warning('Times series a has non-normalizable subsequences, which were not detected by the missing data passes');
end

invn_a(~isvalid_a) = NaN;
invn_b(~isvalid_b) = NaN;

df = [0; (1/2)*(a(1 + w : end) - a(1 : end - w))];
dg = [0; (a(1 + w : end) - mu_a(2 : end)) + (a(1 : end - w) - mu_a(1 : end - 1))];
dx = [0; (1/2)*(b(1 + w : end) - b(1 : end - w))];
dy = [0; (b(1 + w : end) - mu_b(2 : end)) + (b(1 : end - w) - mu_b(1 : end - 1))];

mp_a = repmat(-1, length(a) - w + 1, 1);
mp_b = repmat(-1, length(b) - w + 1, 1);
mp_a(~isvalid_a) = NaN;
mp_b(~isvalid_b) = NaN;
mpi_a = NaN(length(a) - w + 1, 1);
mpi_b = NaN(length(b) - w + 1, 1);

[mp_a, mp_b, mpi_a, mpi_b] = run_mpx(a, b, mu_a, mu_b, invn_a, invn_b, df, dg, dx, dy, mp_a, mp_b, mpi_a, mpi_b, w);
[mp_b, mp_a, mpi_b, mpi_a] = run_mpx(b, a, mu_b, mu_a, invn_b, invn_a, dx, dy, df, dg, mp_b, mp_a, mpi_b, mpi_a, w);


mp_a = sqrt(2 * w * max(0, 1 - mp_a, 'includenan'));
mp_b = sqrt(2 * w * max(0, 1 - mp_b, 'includenan'));

if first_a ~= 1
    mp_a = [NaN(first_a - 1, 1); mp_a];
    mpi_a = [NaN(first_a - 1, 1); mpi_a];
    % offset indices
    mpi_b = mpi_b + (first_a - 1);
end

if first_b ~= 1
    mp_b = [NaN(first_b - 1, 1); mp_b];
    mpi_b = [NaN(first_b - 1, 1); mpi_b];
    % offset indices
    mpi_a = mpi_a + (first_b - 1);
end

if last_a ~= subcount_a
    mp_a = [mp_a; NaN(subcount_a - length(mp_a), 1)]; 
    mpi_a = [mpi_a; NaN(subcount_a - length(mp_a), 1)];
end

if last_b ~= subcount_b
    mp_b = [mp_b; NaN(subcount_b - length(mp_b), 1)]; 
    mpi_b = [mpi_b; NaN(subcount_b - length(mp_b), 1)];
end


end


function [timeseries, isvalidwindow, first, last] = find_valid_windows(timeseries, subseqlen)

% We're going to check for non-normalizable windows here.

% Check for windows containing non-finite data.
subcount = length(timeseries) - subseqlen + 1;
isfinitedata = isfinite(timeseries);
isvalidwindow = movsum(isfinitedata, [0 subseqlen-1], 'Endpoints', 'discard') == subseqlen;

% Typically these come in streams, so we zero them rather than use the
% windowed mean

timeseries(~isfinitedata) = 0;

% Now we need to check for sequences of constants, since these represent
% singularities, which greatly impact stability

issingularity = false(subcount, 1);
% find points with 2 identical consecutive points
for i = 1 : subcount
    % it doesn't matter whether any of these are attributable
    % to zeroing, as this is just a conditioning issue
    issingularity(i) = all(timeseries(i + 1 : i + subseqlen - 1) == timeseries(i));
end

isvalidwindow(issingularity) = false;

% Trim cases where we have bad leading or trailing windows
first = find(isvalidwindow, 1);
last = find(isvalidwindow, 1, 'last');
timeseries = timeseries(first : last + subseqlen - 1);
isvalidwindow = isvalidwindow(first : last);

if isempty(timeseries) || last == first
    return;
end

issingularity = issingularity(first : last);
singularities = find(issingularity);

i = 1;
while i < length(singularities)
    j = i;
    while j < length(singularities)
        % finding the end
        if singularities(j + 1) ~= singularities(j) + 1
            break;
        end
        j = j + 1;
    end
    % find an appropriate scale factor based on the constant occupying
    % this region.
    v = abs(timeseries(singularities(i)));
    if v == 0 
        %just use a standard deviation of 1.
        % It might seem better to look at surrounding windows, but at 
        % that point, we have to recursively check those windows.
        scale = 1;
    elseif v < 1
        % 0 < v < 1, so natural log is negative
        p = 1 + log(v);
        if p == 0
            % 1/2^52 on most systems
            scale = eps(1);
        else
            scale = 1 / abs(p);
        end
    else
        scale = 1 + log(v);
    end
   
    c = 1/64;
    % c is just to ensure this doesn't perturb leading digits
    % in case this data provides a partial window for some normalizable
    % section. Feel free to make it smaller or larger. 
    % The data variance dependent portion is handled by "scale".
    % 1/64 just happens to be a change in exponent only for any number
    % that doesn't hit the de-normal range. I used it for something close
    % to .01 to avoid perturbing leading digits for the most part.
    timeseries(i : j + subseqlen - 1) = timeseries(i : j + subseqlen - 1) + c * scale * randn(j + subseqlen - i, 1);
    
    % If any kind of error happens here, this is the point where we give
    % up, as it becomes very difficult to achieve predictable results.
    if any(~isfinite(timeseries(i : j + subseqlen - 1)))
        error(fprintf('unable to fix missing spanning %d to %d\n', i, j));
    end
    
    i = j + 1;
    
end

end


function [mpa, mpb, mpia, mpib] = run_mpx(a, b, mua, mub, invna, invnb, df, dg, dx, dy, mpa, mpb, mpia, mpib, w)
amx = length(a) - w + 1;
bmx = length(b) - w + 1;

for ia = 1 : amx
    mx = min(amx - ia + 1, bmx);
    c = sum((a(ia : ia + w - 1) - mua(ia)) .* (b(1 : w) - mub(1)));
    for ib = 1 : mx
        c = c + df(ib + ia - 1) * dy(ib) + dg(ib + ia - 1) * dx(ib);
        c_cmp = c * invna(ib + ia - 1) * invnb(ib);
        if c_cmp > mpa(ib + ia - 1)
            mpa(ib + ia - 1) = c_cmp;
            mpia(ib + ia - 1) = ib;
        end
        if c_cmp > mpb(ib)
            mpb(ib) = c_cmp;
            mpib(ib) = ib + ia - 1;
        end
    end
end
end


function [ res ] = moving_mean(a,w)
% moving mean over sequence a with window length w
% based on Ogita et. al, Accurate Sum and Dot Product

% A major source of rounding error is accumulated error in the mean values, so we use this to compensate. 
% While the error bound is still a function of the conditioning of a very long dot product, we have observed 
% a reduction of 3 - 4 digits lost to numerical roundoff when compared to older solutions.

res = zeros(length(a) - w + 1, 1);
p = a(1);
s = 0;

for i = 2 : w
    x = p + a(i);
    z = x - p;
    s = s + ((p - (x - z)) + (a(i) - z));
    p = x;
end

res(1) = p + s;

for i = w + 1 : length(a)
    x = p - a(i - w);
    z = x - p;
    s = s + ((p - (x - z)) - (a(i - w) + z));
    p = x;
    
    x = p + a(i);
    z = x - p;
    s = s + ((p - (x - z)) + (a(i) - z));
    p = x;
    
    res(i - w + 1) = p + s;
end

res = res ./ w;

end

function [dist] = MASS_V2(x, y)
%x is the data, y is the query
m = length(y);
n = length(x);

%compute y stats -- O(n)
meany = mean(y);
sigmay = std(y,1);

%compute x stats -- O(n)
meanx = movmean(x,[m-1 0]);
sigmax = movstd(x,[m-1 0],1);

y = y(end:-1:1);%Reverse the query
y(m+1:n) = 0; %aappend zeros

%The main trick of getting dot products in O(n log n) time
X = fft(x);
Y = fft(y);
Z = X.*Y;
z = ifft(Z);

dist = 2*(m-(z(m:n)-m*meanx(m:n)*meany)./(sigmax(m:n)*sigmay));
dist = sqrt(dist);
dist = dist./sqrt(length(y));

%added 2021-09-01 by Ryan to match standard mp euclidean distances
dist = dist*sqrt(length(x));
end

function [rMatrixProfile, lMatrixProfile, rMatrixProfileIdx, lMatrixProfileIdx] = mpxLeftRight(timeSeries, minlag, subseqLen)



% Code and skew symmetric update formulas by Kaveh Kamgar. They originated
% as a modification to earlier work by Yan Zhu.
% GUI and top k motif critera are set up to match the prior work by Michael Yeh as close as possible.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product
%

n = length(timeSeries);

% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if nargin ~= 3
    error('incorrect number of input arguments');
elseif ~isvector(timeSeries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqLen) && floor(subseqLen) == subseqLen) || (subseqLen < 2) || (subseqLen > length(timeSeries))
    error('subsequence length must be an integer value between 2 and the length of the timeSeries');
end

transposed_ = isrow(timeSeries);
if transposed_
    timeSeries = transpose(timeSeries);
end

nanIdx = isnan(timeSeries);
nanSubseqs = movsum(~isfinite(timeSeries), [0 subseqLen-1], 'Endpoints', 'discard') ~= 0;
rMatrixProfile = repmat(-1, n - subseqLen + 1, 1);
rMatrixProfile(nanSubseqs) = NaN;
lMatrixProfile = repmat(-1, n - subseqLen + 1, 1);
lMatrixProfile(nanSubseqs) = NaN;
timeSeries(nanIdx) = 0;
mu = moving_mean(timeSeries, subseqLen);
invnorm = fastinvn(timeSeries, mu, subseqLen);
invnorm(nanSubseqs) = NaN;


% The update equations diff_f and diff_g are slightly difficult to derive directly, but they arise naturally taking the difference
% between cov(x(i + 1 : i + subseqLen), x(j + 1 : j + subseqLen)) and
% cov(x(i : i + subseqLen - 1), x(j : j + subseqLen - 1))
% from Pebay et al, then applying the identity
% ac - bd = (1/2) * ((a + b) * (b - d) + (a - b) * (b + d))
% to the resulting components.

% Covariance updates which involve a high level of numerical cancellation
% then correspond to very very small changes in covariance and seem to
% rarely ever become a problem in practice, particularly as motif discovery isn't as
% sensitive to perturbations as the algorithms used here to compute cross correlation and euclidean distance.
df = [0; (1/2)*(timeSeries(1 + subseqLen : n) - timeSeries(1 : n - subseqLen))];
dg = [0; (timeSeries(1 + subseqLen : n) - mu(2 : n - subseqLen + 1)) + (timeSeries(1 : n - subseqLen) - mu(1 : n - subseqLen))];

lMatrixProfileIdx = NaN(n - subseqLen + 1, 1);
rMatrixProfileIdx = NaN(n - subseqLen + 1, 1);
comparesTotal = (n - subseqLen + 1) * (n - subseqLen - minlag + 1) / 2;  % <-- each update actually updates 2
updateRate = comparesTotal/100;
comparesSoFar = 0;
counter = 0;

% The terms row and diagonal here refer to a hankel matrix representation of a time series
% This uses scaled cross correlation (scaled to allow the use of movstd) as an intermediate quantity for performance reasons.
% It is later reduced to z-normalized euclidean distance.
for diag = minlag + 1 : n - subseqLen + 1
    if comparesTotal > 32768 && counter >= updateRate
        comparesSoFar = comparesSoFar + counter;
%         fprintf('we are approximately %d percent complete\n', floor(100 * comparesSoFar/comparesTotal));
        counter = 0;
    end
    cov_ = (sum((timeSeries(diag : diag + subseqLen - 1) - mu(diag)) .* (timeSeries(1 : subseqLen) - mu(1))));
    for row = 1 : n - subseqLen - diag + 2
        cov_ = cov_ + df(row) * dg(row + diag - 1) + df(row + diag - 1) * dg(row);
        corr_ = cov_ * invnorm(row) * invnorm(row + diag - 1);
        if corr_ > rMatrixProfile(row)
            rMatrixProfile(row) = corr_;
            rMatrixProfileIdx(row) = row + diag - 1;
        end
        if corr_ > lMatrixProfile(row + diag - 1)
            lMatrixProfile(row + diag - 1) = corr_;
            lMatrixProfileIdx(row + diag - 1) = row;
        end
    end
    counter = counter + (n - subseqLen - diag + 2);
    
end

% This is the correct ordering. findMotifsDiscords uses a correlation based
% updates to avoid it being problematically slow if called in a loop.
%[motifsIdx, mpAugmented] = findMotifs(timeSeries, mu, invnorm, matrixProfile, matrixProfileIdx, subseqLen, 3, 10, minlag, 2);
%[discordsIdx] = findDiscords(mpAugmented, 3, minlag);
%timeSeries(nanIdx) = NaN;
rMatrixProfile = sqrt(2 * subseqLen * (1 - min(1, rMatrixProfile, 'includenan')));
lMatrixProfile = sqrt(2 * subseqLen * (1 - min(1, lMatrixProfile, 'includenan')));


if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    rMatrixProfile = transpose(rMatrixProfile);
    rMatrixProfileIdx = transpose(rMatrixProfileIdx);
    lMatrixProfile = transpose(lMatrixProfile);
    lMatrixProfileIdx = transpose(lMatrixProfileIdx);
end

end

function [invn] = fastinvn(ts, mu, sublen)
% This is a simple variation on Welford's method.
% This version still results in some minor cancellation and could be
% improved. It isn't prone to anything disastrous.
invn = zeros(length(mu), 1);
invn(1) = sum((ts(1 : sublen) - mu(1)).^2);
for i = 2:length(invn)
    invn(i) = invn(i - 1) + ((ts(i - 1) - mu(i - 1)) + (ts(i + sublen - 1) - mu(i))) * (ts(i + sublen - 1) - ts(i - 1));
end
invn = 1./sqrt(invn);
end


function [plato, plato_twin, CP] = ContrastProfile(positiveTS, negativeTS, m, forcePlot)
    %%% This is the "classic" Contrast Profile
    %%% Input:
    %%%   positiveTS: A time series containing at least two instances of a desired behavior
    %%%   negativeTS: A time series containing zero instances of a desired behavior
    %%%   m: The approximate window size of the desired behavior
    %%%   forcePlot: (optional) Display the following two plots
    %%%
    %%% Output:
    %%%   plato: The subsequence that most distinguishes 
    %%%          positiveTS from negativeTS
    %%%   plato_twin: The nearest neighbor of plato within positiveTS
    %%%   CP: Contrast Profile, which indicates subsequence within
    %%%       positiveTS that are less conserved in negativeTS

    initial__pos_subcount = length(positiveTS) - m + 1;
    initial__neg_subcount = length(negativeTS) - m + 1;
    
    if nargin == 3
        forcePlot = false;
    elseif nargin ~= 4
        error('incorrect number of input arguments');
    elseif ~isvector(positiveTS) || ~isvector(negativeTS)
        error('first and second arguments must be a 1D vector');
    elseif ~(isfinite(m) && floor(m) == m) || (m < 2) || (initial__pos_subcount < 2) || (initial__neg_subcount < 2)
        error('subsequence length must be an integer value between 2 and the length of the timeseries');
    end
    
    
    dataOrientation = 0; %0 for row, 1 for column
    %Change to row vector for internal use
    if size(positiveTS,1) == 1 
        %save the data orientation for matching output format to input
        %give orientation priority to positiveTS
        dataOrientation = 1; 
        positiveTS = positiveTS';
    end
    if size(negativeTS,1) == 1 
        %give orientation priority to positiveTS, do not save negative
        %orientation if different than positive
        negativeTS = negativeTS';
    end
    
    
    %%% Matrix profile self-join using positiveTS
    [MP_AA, MP_AA_Indices] = mpxSelfJoin(positiveTS,ceil(m/2),m);
    MP_AA = real(MP_AA);
    %%% Euclidean distance values above sqrt(2*m) are equivalent to
    %%%   anti-correlated values
    MP_AA = clipMatrixProfileAmplitude(MP_AA, m);
    %%% pad with NaN to make future comparisons between matrix profiles
    padLength = length(positiveTS) - length(MP_AA) + 1;
    MP_AA = [MP_AA;NaN(padLength,1)];

    %%% Matrix profile AB-join between positiveTS and negativeTS
    [MP_AB, MP_AB_Indices] = mpx_ABBA_v2(positiveTS, negativeTS, m);
    MP_AB = real(MP_AB);
    %%% Euclidean distance values above sqrt(2*m) are equivalent to
    %%%   anti-correlated values
    MP_AB = clipMatrixProfileAmplitude(MP_AB, m);
    %%% pad with NaN to make future comparisons between matrix profiles
    padLength = length(positiveTS) - length(MP_AB) + 1;
    MP_AB = [MP_AB;NaN(padLength,1)];
    
    %%% Contrast Profile
    CP = MP_AB - MP_AA;
    %%% Normalize values to the range [0,1]
    CP = normalizeContrastProfileAmplitude(CP, m);
    
    
    %%% plato is the subsequence in positiveTS corresponding to index with
    %%%   largest contrast profile value
    [maxContrastValue, platoIndex] = max(CP);
    plato = positiveTS(platoIndex:platoIndex + m - 1);
    
    %%% plato_twin is the nearest neighbor of plato within positieTS
    %%%   It is not necessarily the second largest contrast profile value
    platoTwinIndex = MP_AA_Indices(platoIndex);
    plato_twin = positiveTS(platoTwinIndex:platoTwinIndex + m - 1);
    
    
    

    %%%%%%%%%%%%%
    %%% PLOTS %%%
    %%%%%%%%%%%%%
    if forcePlot == true
        
        redColor = [0.73,0.05,0];
        greenColor = [0,0.73,0.41]; 
        blueColor = [0,0.29,0.73];
        grayColor = [0.75,0.75,0.75];
        lightBlueColor = [0.01, 0.83,0.99];
        platoColor = [129/255, 51/255, 144/255];
        platoTwinColor = [115/255, 170/255, 43/255];

        tsLength = length(positiveTS);
        maxTSLength = max(length(positiveTS),length(negativeTS));

        fig = figure('Name','Contrast Profile: Panel','NumberTitle','off'); 
        set(gcf, 'Position', [0,100,2000,600]);
        
        tiledlayout(5,1);
        
        ax1 = nexttile;
        plot(negativeTS,'Color',redColor);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(-)}\\color{black}: Negative Time Series",redColor(1), redColor(2), redColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(negativeTS)],'ytick',[], 'TickDir','out');
        xlim([1,maxTSLength]);
        box off;
        
        ax2 = nexttile;
        plot(positiveTS);
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(+)}\\color{black}: Positive Time Series",blueColor(1),blueColor(2),blueColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[], 'TickDir','out');
        xlim([1,maxTSLength]);
        box off;
        
        ax3 = nexttile;
        plot(MP_AA,'Color',blueColor);
        hold on;
        plot(MP_AB,'Color',redColor);
        xlim([1,maxTSLength]);
        ylim([0,sqrt(2*m)]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}MP^{(+ -)} AB-join, \\color[rgb]{%f,%f,%f}MP^{(+ +)} Self-Join",redColor(1), redColor(2), redColor(3),blueColor(1),blueColor(2),blueColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,sqrt(2*m)], 'TickDir','out');
        box off;
        
        ax4 = nexttile;
        plot(CP,'Color',grayColor);
        hold on;
        scatter(platoIndex, 1,10,'MarkerFaceColor',platoColor,'MarkerEdgeColor',platoColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
        scatter(platoTwinIndex, 1,10,'MarkerFaceColor',platoTwinColor,'MarkerEdgeColor',platoTwinColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0);
        hold off;
        xlim([1,maxTSLength]);
        ylim([0,1]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Contrast Profile, \\color[rgb]{%f,%f,%f}Plato, \\color[rgb]{%f,%f,%f}Plato Twin",grayColor(1), grayColor(2), grayColor(3), platoColor(1), platoColor(2), platoColor(3), platoTwinColor(1), platoTwinColor(2), platoTwinColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,1], 'TickDir','out');
        box off;
        
        ax5 = nexttile;
        distanceProfile = mpx_ABBA_v2(positiveTS, plato, m);
        distanceProfile = real(distanceProfile);
        plot(distanceProfile);
        xlim([1,maxTSLength]);
        formattedTitle = sprintf("Distance Profile: T^{(+)} join Plato");
        title(formattedTitle);
        set(gca,'xtick',[1,length(positiveTS)],'ytick',[0,sqrt(2*m)], 'TickDir','out');
        box off;
        
        linkaxes([ax1 ax2 ax3 ax4 ax5],'x')
        
        %%%%%%%%%%%%%%%%%%
        %%% PLATO plot %%%
        %%%%%%%%%%%%%%%%%%
        negDistanceProfile = MASS_V2(negativeTS, plato);
        [~,negNNIndex] = min(negDistanceProfile);
        negNN = negativeTS(negNNIndex:negNNIndex+m-1);

        fig = figure('Name','Contrast Profile: Plato','NumberTitle','off');
        set(gcf, 'Position', [0,100,400,400]);
        annotation('textbox', [0, 0.85, 0, 0], 'string', 'Plato');
        annotation('textbox', [0, 0.58, 0, 0], 'string', 'Plato Twin');
        annotation('textbox', [0, 0.3, 0, 0], 'string', 'Neg NN');
        
        subsequences = nan(3,m);
        subsequences(1,:) = positiveTS(platoIndex:platoIndex+m-1);
        subsequences(2,:) = positiveTS(platoTwinIndex:platoTwinIndex+m-1);
        subsequences(3,:) = negNN;
        colors = [platoColor; platoTwinColor; redColor];
        plotIndex = 0;
        inset = 0.9;
%         plot(0,0);
        hold on;
        for i = 1:size(subsequences,1)
            color = colors(i,:);
            tempTS = subsequences(i,:);
            tempMin = min(tempTS);
            tempMax = max(tempTS);
            tempRange = max(1e-5, tempMax-tempMin);
            plot(-plotIndex+inset*(tempTS-tempMin)/tempRange,'Color',color);
            plotIndex = plotIndex + 1;
        end
        hold off;
        

        xlim([0,m]);
        formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Plato \\color{black}(top) and \\color[rgb]{%f,%f,%f}Plato Twin \\color{black}(middle)",platoColor(1), platoColor(2), platoColor(3), platoTwinColor(1), platoTwinColor(2), platoTwinColor(3));
        title(formattedTitle);
        set(gca,'xtick',[1,m],'ytick',[], 'TickDir','out');
        box off;
        
        fprintf('plato index: %d, plato twin index: %d\n', platoIndex, platoTwinIndex);
    end
    
    if dataOrientation == 1
       plato = plato'; 
       plato_twin = plato_twin';
       CP = CP';
    end
    
end


function [matrixProfile, matrixProfileIdx, isvalidwindow, motifsIdx, discordsIdx] = mpxSelfJoin(timeseries, minlag, subseqlen, plot_output)
% (mpx_v3 titled from Kaveh Kamgar, 2020-12-24 12:58)
% matrixProfile - distance between each subsequence timeseries(i : i + subseqLen - 1)
%                 for i = 1 .... length(timeseries) - subseqLen + 1, with the condition
%                 that the nearest neighbor of the ith subsequence is at least minlag steps
%                 away from it.
%
% matrixProfileIdx - If matrixProfileIdx(i) == j, then |i - j| >= minlag
%                    and MatrixProfile(i) == norm(zscore(subseq(i), 1) - norm(zscore(subseq(j), 1));
%
% motifsIdx - A 10 x 3 matrix containing the indices of the top motifs and
%             their neighbors. The first two entries in each column
%             indicate a pair. The remaining entries in that column are their
%             neighbors. This implementation picks motifs from best to
%             worst. It picks neighbors using a greedy strategy with an
%             excluded region applied after each pick to all subsequent
%             choices of motifs and neighbors. Outliers are picked first,
%             followed by the top motif and its neighbors, followed by the
%             second and third motifs in the same manner. This avoids
%             picking outliers as motifs on short time series.
%
%             The distance between a motif and its neighbor must be less than
%             two times the distance between a motif and its nearest
%             neighbor. This is not strict enough to guarantee any causal
%             relation, and in cases of weak matches, neighbors may vary
%             independently from the first pair (that is they may be
%             uncorrelated).
%
% discordsIdx - A 3 x 1 colum indicating the furthest outliers by maximal
%               normalized euclidean distance. These are chosen prior to
%               motifs in the current implementation.
%
%

% Edited 2/7/20: Changed update formula again due to numerical cancellation
% issues in some extreme cases, particularly ill conditioned time series
% with missing data. Indexing feature is still not robust to very poorly
% conditioned data.
%
% Edited: 12/24/20: Added missing data features. These are still somewhat
% experimental. They check for any non-normalizable windows, mark them, 
% and apply noise to any sequences of constants. The extra check of update 
% or no update is implicit, since the comparison is false if either operand
% is nan. This is not advisable in all environments.
%
% The difference equations formula was also adjusted to omit the leading 
% zero, because this casues problems in tiled implementations when starting
% from some point other than the beginning of the time series.

% Code and update formulas are by Kaveh Kamgar.
% GUI and top k motif critera are based on but not identical to some code
% by Michael Yeh. Implementation details are specified above.
%
% Additional References
% Yan Zhu, et al, Matrix Profile II: Exploiting a Novel Algorithm and GPUs to break the one Hundred Million Barrier for Time Series Motifs and Join
% Zachary Zimmerman, et al, Scaling Time Series Motif Discovery with GPUs: Breaking the Quintillion Pairwise Comparisons a Day Barrier. (pending review)
% Philippe Pebay, et al, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Takeshi Ogita, et al, Accurate Sum and Dot Product

initial_subcount = length(timeseries) - subseqlen + 1;

% difference equations have 0 as their first entry here to simplify index
% calculations slightly. Alternatively, it's also possible to swap this to the last element
% and reorder the comparison step (or omit on the last step). This is a
% special case when comparing a single time series to itself. The more general
% case with time series A,B can be computed using difference equations for
% each time series.

if nargin == 3
    plot_output = true;
elseif nargin ~= 4
    error('incorrect number of input arguments');
elseif ~isvector(timeseries)
    error('first argument must be a 1D vector');
elseif ~(isfinite(subseqlen) && floor(subseqlen) == subseqlen) || (subseqlen < 2) || (initial_subcount < 2)
    error('subsequence length must be an integer value between 2 and the length of the timeseries');
end

transposed_ = isrow(timeseries);
if transposed_
    timeseries = transpose(timeseries);
end

[ts, isvalidwindow, first, last] = find_valid_windows_SelfJoin(timeseries, subseqlen);
subcount = length(ts) - subseqlen + 1;
if subcount < minlag
    warning("no valid comparisons could be performed");
    matrixProfile = inf(subcount, 1);
    matrixProfileIdx = repmat(-1, subcount, 1);
    motifsIdx = [];
    discordsIdx = [];
    return;
end

mu = moving_mean(ts, subseqlen);
mus = moving_mean(ts(2:end-1), subseqlen - 1);
invnorm = zeros(subcount, 1);


% Methods fail more often here, including the builtin matlab movstd.
% This uses a simple approach, because it's better at detecting constant
% regions, which should be completely ignored.
for i = 1 : subcount
    invnorm(i) = 1 ./ norm(ts(i : i + subseqlen - 1) - mu(i));
end

if any(~isfinite(invnorm) & isvalidwindow)
    % This is an indicator that all of the preconditioners failed. It may
    % or may not impact the result.
    warning('non finite values encountered when computing reciprocal per window norm.');
end

invnorm(~isvalidwindow) = NaN;
invnorm(~isfinite(invnorm)) = NaN;

% This method requires 4 arrays for rank 1 updates. The previous method used
% 2 arrays in the case of iterating over a single time series. This seems
% to have better stability when dealing with time series with missing data.
% It still loses a number of digits of precision in very badly conditioned
% data.

% References
% Phillippe Pebay, Formulas for Robust, One-Pass Parallel Computation of Covariances and Arbitrary-Order Statistical Moments
% Knuth, semi-numerical algorithms, vol 2

% This was changed to remove the leading zero, as the leading zero
% makes for a bad reference. It doesn't work if you need to tile the
% problem and start from row > 1 rather than row == 1.
dr_bwd = ts(1 : subcount - 1) - mu(1 : subcount - 1);
dc_bwd = ts(1 : subcount - 1) - mus;
dr_fwd = ts(subseqlen + 1 : end) - mu(2 : subcount);
dc_fwd = ts(subseqlen + 1 : end) - mus;
matrixProfile = repmat(-1, subcount, 1);
matrixProfile(~isvalidwindow) = NaN;
matrixProfileIdx = NaN(subcount, 1);

for diag = minlag + 1 : subcount
    cov_ = sum((ts(diag : diag + subseqlen - 1) - mu(diag)) .* (ts(1 : subseqlen) - mu(1)));
    for row = 1 : subcount - diag + 1
        col = diag + row - 1;
        if row > 1
            cov_ = cov_ - dr_bwd(row-1) * dc_bwd(col-1) + dr_fwd(row-1) * dc_fwd(col-1);
        end
        corr_ = cov_ * invnorm(row) * invnorm(col);
        if corr_ > matrixProfile(row)
            matrixProfile(row) = corr_;
            matrixProfileIdx(row) = col;
        end
        if corr_ > matrixProfile(col)
            matrixProfile(col) = corr_;
            matrixProfileIdx(col) = row;
        end
    end
end


% updated to pick outliers independent of motifs
% this means they might overlap or whatever. Compute more of them if you
% need to in order to avoid this. The other methods introduce bias.
[discordsIdx] = findDiscords(matrixProfile, minlag);

[motifsIdx] = findMotifs(ts, mu, invnorm, matrixProfile, matrixProfileIdx, subseqlen, minlag);


% Max caps anything that rounds inappropriately. This can hide potential
% issues at times, but it's fairly rare.
matrixProfile = sqrt(max(0, 2 * subseqlen * (1 - matrixProfile), 'includenan'));

% expand initial
if first ~= 1
    isvalidwindow = [false(first - 1, 1); isvalidwindow];
    matrixProfile = [repmat(-1, first - 1, 1); matrixProfile];
    % offset index to compensate for the dropped leading windows
    matrixProfileIdx = [repmat(-1, first - 1, 1); matrixProfileIdx + first - 1];
end

if last ~= initial_subcount
    extendby = initial_subcount - last + 1;
    isvalidwindow = [isvalidwindow; false(extendby, 1)];
    matrixProfile = [matrixProfile; NaN(extendby, 1)];
    matrixProfileIdx = [matrixProfileIdx; NaN(extendby, 1);];
end

% if plot_output
%     gui = mpgui_ed(timeseries, subseqlen);
%     % first pair is plottted alongside data
%     gui.plotProfile(matrixProfile);
%     if isfinite(motifsIdx(1, 1)) && ~isnan(motifsIdx(2, 1))
%         gui.plotData(motifsIdx(1, 1), motifsIdx(2, 1));
%         motiftitle = sprintf('The first motif pair is located at %d and %d.', motifsIdx(1,1), motifsIdx(2,1));
%         gui.plotMotif(gui.motifAx1, motifsIdx(1, 1), motifsIdx(2:end, 1), motiftitle);
%         if isfinite(motifsIdx(1, 2))
%             motiftitle = sprintf('The second motif pair is located at %d and %d.', motifsIdx(1,2), motifsIdx(2,2));
%             gui.plotMotif(gui.motifAx2, motifsIdx(1, 2), motifsIdx(2:end, 2), motiftitle);
%         end
%         if isfinite(motifsIdx(1, 3))
%             motiftitle = sprintf('The third motif pair is located at %d and %d.', motifsIdx(1,3), motifsIdx(2,3));
%             gui.plotMotif(gui.motifAx3, motifsIdx(1, 3), motifsIdx(2:end, 3), motiftitle);
%         end
%         gui.plotDiscords(discordsIdx, sprintf(['The top three discords ', '%d(blue), %d(red), %d(green)'], discordsIdx(1), discordsIdx(2), discordsIdx(3)));
%     else
%         gui.plotData();
%     end
%     
%     gui.drawgui;
%     
% end

if transposed_  % matches the profile and profile index but not the motif or discord index to the input format
    matrixProfile = transpose(matrixProfile);
    matrixProfileIdx = transpose(matrixProfileIdx);
end

end

function [timeseries, isvalidwindow, first, last] = find_valid_windows_SelfJoin(timeseries, subseqlen)

% We're going to check for non-normalizable windows here.

% Check for windows containing non-finite data.
subcount = length(timeseries) - subseqlen + 1;
isfinitedata = isfinite(timeseries);
isvalidwindow = movsum(isfinitedata, [0 subseqlen-1], 'Endpoints', 'discard') == subseqlen;

% Typically these come in streams, so we zero them rather than use the
% windowed mean

timeseries(~isfinitedata) = 0;

% Now we need to check for sequences of constants, since these represent
% singularities, which greatly impact stability

issingularity = false(subcount, 1);
% find points with 2 identical consecutive points
for i = 1 : subcount
    % it doesn't matter whether any of these are attributable
    % to zeroing, as this is just a conditioning issue
    issingularity(i) = all(timeseries(i + 1 : i + subseqlen - 1) == timeseries(i));
end

isvalidwindow(issingularity) = false;

% Trim cases where we have bad leading or trailing windows
first = find(isvalidwindow, 1);
last = find(isvalidwindow, 1, 'last');
timeseries = timeseries(first : last + subseqlen - 1);
isvalidwindow = isvalidwindow(first : last);

if isempty(timeseries) || last == first
    return;
end

issingularity = issingularity(first : last);
singularities = find(issingularity);

i = 1;
while i < length(singularities)
    j = i;
    while j < length(singularities)
        % finding the end
        if singularities(j + 1) ~= singularities(j) + 1
            break;
        end
        j = j + 1;
    end
    % find an appropriate scale factor based on the constant occupying
    % this region.
    v = abs(timeseries(singularities(i)));
    if v == 0 
        %just use a standard deviation of 1.
        % It might seem better to look at surrounding windows, but at 
        % that point, we have to recursively check those windows.
        scale = 1;
    elseif v < 1
        % 0 < v < 1, so natural log is negative
        p = 1 + log(v);
        if p == 0
            % 1/2^52 on most systems
            scale = eps(1);
        else
            scale = 1 / abs(p);
        end
    else
        scale = 1 + log(v);
    end
   
    % 0.01 is just to ensure this doesn't perturb leading digits
    % in case this data provides a partial window for some normalizable
    % section. Feel free to make it smaller or larger. 
    % The data variance dependent portion is handled by "scale".
    timeseries(i : j + subseqlen - 1) = 0.01 * scale * randn(j + subseqlen - i, 1);
    
    % If any kind of error happens here, this is the point where we give
    % up, as it becomes very difficult to achieve predictable results.
    if any(~isfinite(timeseries(i : j + subseqlen - 1)))
        error(fprintf('unable to fix missing spanning %d to %d\n', i, j));
    end
    
    i = j + 1;
    
end

end


function [discordIdx] = findDiscords(matrixProfile, exclusionLen)
% This function assumes a correlation based implementation of matrixProfile
% and that any excluded regions have been set to NaN, nothing else.

% Changed on 1/29/20
% Making motif discovery and discord discovery co-dependent introduces odd
% issues. I am avoiding any sharing of excluded regions at this point,
% since it's not clear which should come first. Computing 3 or k or
% whatever motifs prior to discords was definitely not correct, but I'm not
% sure this is a better idea.

discordCount = 3;
discordIdx = zeros(discordCount, 1);
[~, idx] = sort(matrixProfile);

for i = 1 : discordCount
    f = find(isfinite(matrixProfile(idx)) & (matrixProfile(idx) > -1), 1);
    if isempty(f) || (matrixProfile(f) == -1)
        discordIdx(i : end) = NaN;
        break;
    end
    discordIdx(i) = idx(f);
    exclRangeBegin = max(1, discordIdx(i) - exclusionLen + 1);
    exclRangeEnd = min(length(matrixProfile), discordIdx(i) + exclusionLen - 1);
    matrixProfile(exclRangeBegin : exclRangeEnd) = NaN;
end
end


function [motifIdxs, matrixProfile] = findMotifs(timeseries, mu, invnorm, matrixProfile, profileIndex, subseqLen, exclusionLen)
% This is adapted match the output of some inline code written by Michael Yeh
% to find the top k motifs in a time series. Due to some bug fixes I applied, the two
% may return slightly different results, although they almost always agree on the top 2 motif pairs.


% +2 allows us to pack the original motif pair with its neighbors.
% Annoyingly if it's extended, matlab will pad with zeros, leading to some
% rather interesting issues later.
motifCount = 3;
radius = 2;
neighborCount = 10;

motifIdxs = NaN(neighborCount + 2, motifCount);

% removing use of ffts since this formula is easier to analyze and conv2
% is pretty well optimized
crosscov = @(idx) conv2(timeseries, (timeseries(idx + subseqLen - 1 : -1 : idx) - mu(idx)) .* invnorm(idx), 'valid');

for i = 1 : motifCount
    [corr_, motIdx] = max(matrixProfile);
    % -1 means this is maximally negatively correlated and was therefore
    % never updated
    if ~isfinite(corr_) || corr_ == -1
        break;
    end
    % order subsequence motif pair as [time series index of 1st appearance, time series index of 2nd appearance]
    motifIdxs(1 : 2, i) = [min(motIdx, profileIndex(motIdx)), max(motIdx, profileIndex(motIdx))];
    [corrProfile] = crosscov(motIdx);
    corrProfile = min(1, corrProfile(1 : length(timeseries) - subseqLen + 1) .* invnorm, 'includenan');
    % This uses correlation instead of normalized Euclidean distance, because it's easier to work with and involves fewer operations.
    
    corrProfile(isnan(matrixProfile)) = NaN;
    if exclusionLen > 0
        for j = 1 : 2
            exclRangeBegin = max(1, motifIdxs(j, i) - exclusionLen + 1);
            exclRangeEnd = min(length(matrixProfile), motifIdxs(j, i) + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    
    for j = 3 : neighborCount + 2
        [neighborCorr, neighbor] = max(corrProfile);
        % If you want to put a reasonable global bound on it,
        % set the if statement to also skip anything where neighborCorr <= 0 as well.
        % This was reverted to match earlier code, which did not enforce such a restriction.
        %
        if  ~isfinite(neighborCorr) || ((1 - neighborCorr) >= radius * (1 - corr_))
            break;
        end
        motifIdxs(j, i) = neighbor;
        if exclusionLen > 0
            exclRangeBegin = max(1, neighbor - exclusionLen + 1);
            exclRangeEnd = min(length(matrixProfile), neighbor + exclusionLen - 1);
            corrProfile(exclRangeBegin : exclRangeEnd) = NaN;
        end
    end
    % Matlab allows NaNs to be ignored from min/max calculations by default, so this
    % is a Matlab specific way to iteratively add excluded ranges of elements.
    matrixProfile(isnan(corrProfile)) = NaN;
end
end

