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
%%%   noveletPrimaryContrast: Contrast value of each novelet in the K=1 obj.EP
%%%   noveletNaryContrast: Contrast value of each novelet after appending
%%%     the previous Novelet to obj.negativeTS. May be helpful in identifying
%%%     diminishing returns and retundant behaviors.
classdef OnlineEmergenceProfile
    properties
        positiveTS
        negativeTS
        mm % Subsequence length
        startLocalOffset
        endLocalOffset

        exclusionLength
        contextLength
        threshold
        MP_AA
        MP_AA_Indices
        MP_AB
        DP_Novelets
        EP %Highlights seconds instance of behavior

        minBufferLength
        oldNewsAnnotation
        firstFoundAnnotation
        secondFoundAnnotation
        bufferStartIndex

        novelets
        noveletIndices
        noveletNNIndices
        noveletScores
        noveletDistanceThresholds
    end
    methods
        function obj = OnlineEmergenceProfile(positiveTS, negativeTS, mm, threshold, startLocalOffset, endLocalOffset)
            if nargin < 1
               positiveTS = [];
            end
            if nargin < 2 || length(negativeTS) < mm
                negativeTS = [];
            end
            if nargin < 3
                mm = 10;
            end
            if nargin < 4
                threshold = 0;
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
            obj.threshold = threshold;
            
            obj.minBufferLength = obj.mm*4;
            obj.bufferStartIndex = 1;
            obj.oldNewsAnnotation = false(length(obj.positiveTS),1); %%% boolean
            obj.novelets = []; %%% dimension: (numNovelets, mm)
            obj.noveletIndices = [];
            obj.noveletNNIndices = [];
            obj.noveletScores = [];
            obj.noveletDistanceThresholds = [];

            if isempty(obj.negativeTS)
                obj.negativeTS = rand(obj.mm+1,1);
            end

            
            
            obj.bufferStartIndex = 1;

            obj.positiveTS = [];
            obj.MP_AA = [];
            obj.MP_AA_Indices = [];
            obj.MP_AB = [];
            obj.DP_Novelets = [];
            obj = obj.update(positiveTS);

            end

        function obj = discoverNovelets(obj, t, MP_AB, MP_AA, MP_AA_Indices)

            while obj.threshold <= 0 || obj.threshold >= 1
                tempEP = MP_AB - MP_AA;
                tempEP = normalizeContrastProfileAmplitude(tempEP, obj.mm);
                
                obj.plotEmergenceProfile(tempEP);

                [~, peakValues] = exclusionZonePeaks(tempEP, obj.mm, obj.exclusionLength);
                figure;
                histogram(peakValues,100);
                title("Histogram of Subsequence Novelty Scores")
                userMessage = "";
                userMessage = userMessage + "Please choose a novelty threshold within (0,1).\nLarger values require higher novelty.\n";
                obj.threshold = input(userMessage);
            end

            t = reshape(t, length(t), 1);

            
            bufferLength = length(t);
            initialBufferStartIndex = obj.bufferStartIndex;
            while bufferLength >= obj.minBufferLength
                bufferIndex = obj.bufferStartIndex - initialBufferStartIndex + 1;
                startIndexActionWindow = bufferIndex;
                endIndexActionWindow = min(length(obj.positiveTS), startIndexActionWindow + obj.minBufferLength - 1);
                endIndexMPActionWindow = endIndexActionWindow - obj.mm + 1;
                tSub = t(startIndexActionWindow:endIndexActionWindow);
                MP_ABSub = MP_AB(startIndexActionWindow:endIndexMPActionWindow);
                MP_ABSubUpdated = MP_ABSub;
                DP_NoveletsUpdated = 2*sqrt(obj.mm)*ones(length(MP_ABSub),1);
                MP_AASub = MP_AA(startIndexActionWindow:endIndexMPActionWindow);
                MP_AA_IndicesSub = MP_AA_Indices(startIndexActionWindow:endIndexMPActionWindow);

                %%% Update the MP_AB here
                
                for jj = 1:size(obj.novelets,1)
                    try
                        novelet = obj.novelets(jj, :)';
                        newDP_Novelets = mpx_ABBA_v2(tSub, novelet, obj.mm); %%%TODO: verify if this overlap is correct
                        newDP_Novelets = clipMatrixProfileAmplitude(newDP_Novelets, obj.mm); 
                        DP_NoveletsUpdated(1:length(newDP_Novelets)) = min(DP_NoveletsUpdated(1:length(newDP_Novelets)), newDP_Novelets,'includenan');
                    catch
%                         fprintf("Error when computing Novelet distance profile\n");
%                         return
                    end
                 end
                
                DP_NoveletsUpdated = min(sqrt(2*obj.mm)*ones(length(DP_NoveletsUpdated), 1), DP_NoveletsUpdated);

                %%%Ignore time series nearby
                if obj.bufferStartIndex > obj.startLocalOffset && obj.startLocalOffset > 0
                    startLocalOffset = obj.bufferStartIndex - obj.startLocalOffset;
                    endLocalOffset = obj.bufferStartIndex - obj.endLocalOffset;
                    localTS = obj.positiveTS(startLocalOffset: endLocalOffset);
                    newDP_Novelets = mpx_ABBA_v2(tSub, localTS, obj.mm); %%%TODO: verify if this overlap is correct
                    newDP_Novelets = clipMatrixProfileAmplitude(newDP_Novelets, obj.mm); 
                    DP_NoveletsUpdated(1:length(newDP_Novelets)) = min(DP_NoveletsUpdated(1:length(newDP_Novelets)), newDP_Novelets);
                end
                

                MP_ABSubUpdated = min(MP_ABSub, DP_NoveletsUpdated, 'includenan');

                EPSub = MP_ABSubUpdated - MP_AASub;
                EPSub = normalizeContrastProfileAmplitude(EPSub, obj.mm);
                if any(EPSub>0)
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
                   [peakIndices, peakValues] = exclusionZonePeaks(EPSub, obj.mm, obj.exclusionLength, 2);                   
                   if length(peakIndices) >= 2
                       peakIndex = peakIndices(1);
                   else
                        peakIndex = 1;
                   end
                else
                    peakIndex = 1;
                end
                
                peakMP_AB = MP_ABSubUpdated(peakIndex);
                peakMP_AA = MP_AASub(peakIndex);
                peakMP_AA_Index = MP_AA_IndicesSub(peakIndex);
% 
%                 
%                 peakContrast = peakMP_AB - peakMP_AA;
%                 peakContrast = normalizeContrastProfileAmplitude(peakContrast, obj.mm);
%                   peakMP_AA_Index = MP_AA_IndicesSub(peakIndex);
                  peakContrast = EPSub(peakIndex);

                if peakContrast > obj.threshold
                    obj.noveletNNIndices(end+1) = obj.bufferStartIndex + peakIndex - 1; %%% Swapped because we find the second instance
                    obj.noveletIndices(end+1) = peakMP_AA_Index;                        %%%   then look back to the first
                    novelet = obj.getSubsequenceWithContext(peakMP_AA_Index, obj.contextLength);
                    obj.novelets(end+1,:) = novelet;
                    
                    obj.noveletScores(end+1) = peakContrast;
%                     obj.noveletDistanceThresholds(end+1) = (peakMP_AB + peakMP_AA)/2;
                    %%% The next novelet will need to exceed the AA
                    %%% distance + the novelty threshold
                    obj.noveletDistanceThresholds(end+1) = peakMP_AA + 0.99*obj.threshold;


                    startIndex = obj.bufferStartIndex + peakIndex - 1;
                    endIndex = startIndex + obj.mm - 1;
                    obj.oldNewsAnnotation(startIndex: endIndex) = true;

                    startIndex = peakMP_AA_Index;
                    endIndex = startIndex + obj.mm - 1;
                    obj.oldNewsAnnotation(startIndex: endIndex) = true;
                end

                
                activeEndIndex = obj.bufferStartIndex + length(MP_AASub)-1;

                obj.EP(obj.bufferStartIndex:activeEndIndex) = EPSub;

                
                obj.bufferStartIndex = obj.bufferStartIndex + peakIndex + obj.exclusionLength;
                bufferLength = length(t) - (obj.bufferStartIndex - initialBufferStartIndex + 1) + 1;
            end
        end

        function obj = update(obj, t)
            if isempty(t)
                return;
            else
                t = reshape(t, length(t), 1);
                obj.positiveTS = [obj.positiveTS; t];
                obj.EP = [obj.EP; zeros(length(t),1)];
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
                tempMP_AB = real(tempMP_AB);
            else
                tempMP_AB = [];
            end

            newMP_AB(1:length(tempMP_AB)) = tempMP_AB;
            %%% Start Novelet Comparisons
            %%%TODO: remove, novelet comparison are done in novelet
            %%%discovery function
%             numNovelets = size(obj.novelets,1);
%             noveletDP = sqrt(2*obj.mm)*ones(numNovelets,unprocessedLength);
%             
%             if numNovelets > 0
%                 for noveletIndex = 1:numNovelets
%                     novelet = reshape(obj.novelets(noveletIndex,:), obj.mm + 2*obj.contextLength, 1);
%                     tempNoveletDist = real(MASS_V2(unprocessedPositiveTS, novelet));
%                     tempNoveletDist = clipMatrixProfileAmplitude(tempNoveletDist, obj.mm);
%                     noveletDP(noveletIndex,1:length(tempNoveletDist)) = tempNoveletDist;
%                 end
%                 [noveletMinDist, ~] = min(noveletDP,[],1);
%                 noveletMinDist = reshape(noveletMinDist, unprocessedLength, 1);
%     
%                 newMP_AB = min(newMP_AB, noveletMinDist);
%             end


            %%%Extending the left only matrix profile can be done by
            %%%combining an AB-join between unprocessed posTS with rest of
            %%%posTS , then the left only self-join of the unprocessed
            %%%posTS
            
            newMP_AA = sqrt(2*obj.mm)*ones(unprocessedLength, 1);
            newMP_AA_Indices = nan(unprocessedLength, 1);

%             [~, tempLeftAA, ~, tempLeftAAIndices] = mpxLeftRight(obj.positiveTS, ceil(obj.mm/2), obj.mm);
%             tempLeftAA = tempLeftAA(obj.bufferStartIndex:end);
%             newMP_AA(1:length(tempLeftAA)) = tempLeftAA;
            tempAB = nan(unprocessedLength,1);
            if obj.bufferStartIndex > obj.mm
                [tempAB, ~, tempABIndices] = mpx_ABBA_v2(unprocessedPositiveTS, obj.positiveTS(1:obj.bufferStartIndex-1), obj.mm);
            end
            tempAB = real(tempAB);
            [~, tempLeftAA, ~, tempLeftAAIndices] = mpxLeftRight(unprocessedPositiveTS, obj.exclusionLength, obj.mm);
            tempLeftAA = real(tempLeftAA);
            tempLeftAAIndices = tempLeftAAIndices + obj.bufferStartIndex -1;
            tempLeftAA(1:2*obj.mm) = sqrt(2*obj.mm)*ones(2*obj.mm,1);
            
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
            
            obj.MP_AA(obj.bufferStartIndex:end) = [];
            obj.MP_AA_Indices(obj.bufferStartIndex:end) = [];
            obj.MP_AB(obj.bufferStartIndex:end) = [];

            obj.MP_AA = [obj.MP_AA; newMP_AA];
            obj.MP_AA_Indices = [obj.MP_AA_Indices; newMP_AA_Indices];
            obj.MP_AB = [obj.MP_AB; newMP_AB];
 

            

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

        function obj = resetThreshold(obj, threshold)
            %%% Set new threshold
            obj.threshold = threshold;
            
            %%% Reset relevant variables
            obj.bufferStartIndex = 1;
            obj.oldNewsAnnotation = false(length(obj.positiveTS),1); %%% boolean
            obj.novelets = []; %%% dimension: (numNovelets, mm)
            obj.noveletIndices = [];
            obj.noveletNNIndices = [];
            obj.noveletScores = [];
            obj.noveletDistanceThresholds = [];

            %%% Rediscover Novelets
            obj = discoverNovelets(obj, obj.positiveTS, obj.MP_AB, obj.MP_AA, obj.MP_AA_Indices);
        end

        function f = plotEmergenceProfile(obj, ep)
            %%%%%%%%%%%%%
            %%% PLOTS %%%
            %%%%%%%%%%%%%    

            redColor = [0.73,0.05,0];
%             greenColor = [0,0.73,0.41]; 
            blueColor = [0,0.29,0.73];
            negColor = [130/255,0,0];
%             greenColor = [0,0.73,0.41]; 
            posColor = [0,32/255,96/255];
            grayColor = [0.75,0.75,0.75];
            lightGrayColor = [0.9, 0.9, 0.9];
%             lightBlueColor = [0.01, 0.83,0.99];
            emergenceProfileColor = [0, 176/255, 80/255];
            noveletColor = [143/255, 226/255, 227/255];
            noveletLabelColor = [56/255, 226/255, 235/255];
            noveletNNColor = [227/255, 143/255, 219/255];
            noveletNNLabelColor = [235/255, 56/255, 217/255];
    
%             tsLength = length(obj.positiveTS);
            maxTSLength = max(length(obj.positiveTS),length(obj.negativeTS));
    
            f = figure('Name','Emergence Profile: Panel','NumberTitle','off'); 
            set(gcf, 'Position', [0,100,2000,600]);
            
            tiledlayout(4,1);
            
            ax1 = nexttile;
            if ~isempty(obj.negativeTS)
                plot(obj.negativeTS,'Color',negColor);
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
            formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}T^{(+)}\\color{black}: Positive Time Series",posColor(1),posColor(2),posColor(3));
            title(formattedTitle);
            set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[], 'TickDir','out');
            xlim([1,maxTSLength]);
            box off;
            
            ax3 = nexttile;
            plot(obj.MP_AA,'Color',blueColor);
            hold on;
            plot(obj.MP_AB,'Color',redColor);
            xlim([1,maxTSLength]);
            ylim([0,1.1*sqrt(2*obj.mm)]);
            formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}MP^{(+ -)} AB-join, \\color[rgb]{%f,%f,%f}LMP^{(+ +)} Left Self-Join",redColor(1), redColor(2), redColor(3),blueColor(1),blueColor(2),blueColor(3));
            title(formattedTitle);
            set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[0,sqrt(2*obj.mm)], 'TickDir','out');
            box off;
            
            ax4 = nexttile;
            hold on;
            plot([0,maxTSLength],[obj.threshold, obj.threshold],'--','Color',grayColor);
            plot(ep,'Color',emergenceProfileColor);
            scatter(obj.noveletIndices, ones(length(obj.noveletScores),1) ,20,'v','MarkerFaceColor',noveletLabelColor,'MarkerEdgeColor',noveletLabelColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
            scatter(obj.noveletNNIndices, ones(length(obj.noveletScores),1), 20,'v','MarkerFaceColor',noveletNNLabelColor,'MarkerEdgeColor',noveletNNLabelColor,'MarkerFaceAlpha',0.7,'MarkerEdgeAlpha',0); 
            hold off;
            xlim([1,maxTSLength]);
            ylim([-0.1,1.1]);
            formattedTitle = sprintf("\\color[rgb]{%f,%f,%f}Emergence Profile, \\color[rgb]{%f,%f,%f}First, \\color[rgb]{%f,%f,%f}Second",emergenceProfileColor(1), emergenceProfileColor(2), emergenceProfileColor(3), noveletLabelColor(1), noveletLabelColor(2), noveletLabelColor(3), noveletNNLabelColor(1), noveletNNLabelColor(2), noveletNNLabelColor(3));
            title(formattedTitle);
            set(gca,'xtick',[1,length(obj.positiveTS)],'ytick',[0,1], 'TickDir','out');
            box off;
            
            
            linkaxes([ax1 ax2 ax3 ax4],'x')
        end

        function f = plotNovelets(obj, sortMode)
            if nargin <= 1
                %%% Options to sort by. First is at the top.
                %%% "descendingScore": plot novelets with highest novety
                %%%     score at the top
                %%% "ascendingChronology": plot novelets in order
                %%%     discovered
                sortMode = "ascendingChronology";
            end

            if sortMode == "descendingScore"
                [~,sortedOrder] = sort(obj.noveletScores,"descend");
                sortedIndices = obj.noveletIndices(sortedOrder);
                sortedNNIndices = obj.noveletNNIndices(sortedOrder);
                sortTitle = "Descending Novelty Score";
            else %%% "ascendingChronology"
                sortedOrder = 1:length(obj.noveletIndices);
                sortedIndices = obj.noveletIndices;
                sortedNNIndices = obj.noveletNNIndices;
                sortTitle = "Order Discovered";
            end

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
                    
            %%%Generating unique colors for classes when there can be many classes
            %%% I will assume no more than 1000 results
            numColors = 1000;
            colors = lines(numColors);
            
            f = figure;
            hold on;
            numNovelets = size(obj.novelets,1);
            startIndex = obj.contextLength+1;
            endIndex = startIndex + obj.mm - 1;
            for ki = 1:numNovelets
                %%%Nearest neighbor plotted under
                startIndexNN = sortedNNIndices(ki);
                endIndexNN = startIndexNN + obj.mm - 1;
                tempTS = obj.positiveTS(startIndexNN: endIndexNN);
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
                plot(1:obj.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color', 1-0.3*(1-colors(ki,:)));

                %%% prepare the Novelet without context
                ni = sortedOrder(ki);
                tempTS = obj.novelets(ni,startIndex:endIndex);
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
    
                tempContextTS = obj.novelets(ni,:);
                %%% plot the context in the middle
                plot(-obj.contextLength+1:obj.mm+obj.contextLength,-ki + 0.9*(tempContextTS - tempMin)/tempRange,'Color',lightGrayColor);
                %%% plot the Novelet on top
                plot(1:obj.mm,-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));
            end
            hold off;
            formattedTitle = sprintf("%d Novelets > %.2f, Sorted by %s", numNovelets, obj.threshold, sortTitle);
            title(formattedTitle);
            set(gca,'xtick',[1,obj.mm],'ytick',[], 'TickDir','out');
            xlim([-obj.contextLength+1,obj.mm+obj.contextLength]);
            box off;
        end

        function f = plotNoveletDistances(obj, pos, sortMode)
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
            
            splitPos = [0.2 0.8];
            emptyPosition = [0 splitPos(2) splitPos(1) 1-splitPos(2)];
            posPosition = [splitPos(1) splitPos(2) 1-splitPos(1) 1-splitPos(2)];
            noveletPosition = [0 0 splitPos(1) splitPos(2)];
            distProfPosition = [splitPos(1) 0 1-splitPos(1) splitPos(2)];
            
            neg = obj.negativeTS;
            if nargin <= 1
                pos = obj.positiveTS;
            end
            if nargin <= 2
                %%% Options to sort by. First is at the top.
                %%% "descendingScore": plot novelets with highest novety
                %%%     score at the top
                %%% "ascendingChronology": plot novelets in order
                %%%     discovered
                sortMode = "ascendingChronology";
            end

            f = figure;
            ax1 = subplot('Position', posPosition);
            plot(pos);
            set(gca,'xtick',[],'ytick',[], 'TickDir','out');
            xlim([1,length(pos)]);
            box off;
            
            
            
            if sortMode == "descendingScore"
                [~,sortedOrder] = sort(obj.noveletScores,"descend");
                sortedIndices = obj.noveletIndices(sortedOrder);
                sortedNNIndices = obj.noveletNNIndices(sortedOrder);
                sortTitle = "Descending Novelty Score";
            else %%% "ascendingChronology"
                sortedOrder = 1:length(obj.noveletIndices);
                sortedIndices = obj.noveletIndices;
                sortedNNIndices = obj.noveletNNIndices;
                sortTitle = "Order Discovered";
            end
            
            
                    
            %%%Generating unique colors for classes when there can be many classes
            %%% I will assume no more than 1000 results
            numColors = 1000;
            colors = lines(numColors);
            
            ax2 = subplot('Position', noveletPosition);
            hold on;
            numNovelets = size(obj.novelets,1);
            startIndex = obj.contextLength+1;
            endIndex = startIndex + obj.mm - 1;
            for ki = 1:numNovelets
                ni = sortedOrder(ki);
                novelet = obj.novelets(ni,startIndex:endIndex)';

                %%% Nearest neighbor in negativeTS
                [~,startIndexNN] = min(real(MASS_V2_nan(neg, novelet)));
                endIndexNN = startIndexNN + obj.mm - 1;
                tempTS = obj.negativeTS(startIndexNN:endIndexNN);
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
                plot(-ki + 0.9*(tempTS - tempMin)/tempRange,'Color', [1, 0.95, 0.95]);

                %%%Nearest neighbor plotted under
                startIndexNN = sortedNNIndices(ki);
                endIndexNN = startIndexNN + obj.mm - 1;
                tempTS = obj.positiveTS(startIndexNN: endIndexNN);
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
                plot(-ki + 0.9*(tempTS - tempMin)/tempRange,'Color', 1-0.5*(1-colors(ki,:)));
            
                %%% prepare the Novelet without context
                tempTS = novelet;
                tempMin = min(tempTS);
                tempMax = max(tempTS);
                tempRange = tempMax-tempMin;
            
                %%% plot the Novelet on top
                plot(-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));
            end
            hold off;
            % formattedTitle = sprintf("%d Novelets > %.2f, Sorted by %s", numNovelets, obj.threshold, sortTitle);
            % title(formattedTitle);
            set(gca,'xtick',[1,obj.mm],'ytick',[], 'TickDir','out');
            xlim([1,obj.mm]);
            box off;
            
            ax3 = subplot('Position', distProfPosition);
            numNovelets = size(obj.novelets,1);
            hold on;
            for ki = 1:numNovelets
                %%% prepare the Novelet without context
                ni = sortedOrder(ki);
                startIndex = sortedIndices(ki);
                endIndex = startIndex + obj.mm - 1;
                novelet = obj.positiveTS(startIndex:endIndex);
                tempTS = real(MASS_V2_nan(pos, novelet));
                tempMin = 0;%min(tempTS);
                tempMax = 2*sqrt(length(novelet));%max(tempTS);
                tempRange = tempMax-tempMin;
            
                %%% plot the distance threshold used to classify the behavior
                threshold = obj.noveletDistanceThresholds(ni);
                offsetThreshold = -ki + 0.9*(threshold-tempMin)/tempRange;
                plot([0,length(pos)],[offsetThreshold, offsetThreshold],'--','Color',[0.7, 0.7, 0.7]);
                %%% plot the Novelet on top
                plot(-ki + 0.9*(tempTS - tempMin)/tempRange,'Color',colors(ki,:));
            
                [peakIndices, peakValues] = exclusionZonePeaks(-tempTS, obj.mm, obj.mm, length(tempTS), -threshold);
                for pi = 1:length(peakIndices)
                    offsetMarker = -ki + 0.9*1;
                    scatter(peakIndices, offsetMarker*ones(length(peakIndices),1),'v','filled','MarkerEdgeColor', redColor,'MarkerFaceColor', redColor);
                end
            
            end
            hold off;
            % formattedTitle = sprintf("%d Novelets > %.2f, Sorted by %s", numNovelets, obj.threshold, sortTitle);
            % title(formattedTitle);
            set(gca,'xtick',[],'ytick',[], 'TickDir','out');
            xlim([1,length(pos)]);
            box off;
            
            linkaxes([ax1, ax3],'x');
            linkaxes([ax2 ax3],'y');
        end

        function plot(obj)
            obj.plotEmergenceProfile(obj.EP);
            obj.plotNovelets();
        end

        function savePlot(obj, outputPath, name)
            fig = obj.plotEmergenceProfile(obj.EP);

            fileName = name + "_EmergenceProfile";
            filePath = fullfile(outputPath, fileName + ".fig");
            savefig(fig, filePath);
    
            filePath = fullfile(outputPath, fileName + ".emf");
            print(filePath,'-dmeta');
            close(fig);


            fig = obj.plotNovelets();

            fileName = name + "_Novelets";
            filePath = fullfile(outputPath, fileName + ".fig");
            savefig(fig, filePath);
    
            filePath = fullfile(outputPath, fileName + ".emf");
            print(filePath,'-dmeta');
            close(fig);
        end

        function noveletsWithoutContext = getNoveletsWithoutContext(obj)
            noveletsWithoutContext = obj.novelets(:,obj.contextLength+1:obj.contextLength+obj.mm);
        end

        function obj = setupDummyData(obj, posLength, negLength, numNovelets)
            obj.positiveTS = randn(posLength, 1);
            obj.negativeTS = randn(negLength, 1);
            obj.bufferStartIndex = posLength-obj.minBufferLength+1;
            obj.novelets = randn(numNovelets, obj.mm);
            obj.noveletIndices = ceil(rand(numNovelets,1)*(posLength-3*obj.mm+1));
            obj.noveletNNIndices = ceil(rand(numNovelets,1)*(posLength-3*obj.mm+1));

            obj.noveletScores = randn(numNovelets, 1);
            obj.noveletDistanceThresholds = randn(numNovelets, 1);

            obj.MP_AA = randn(posLength, 1);
            obj.MP_AA_Indices = (1:posLength)';
            obj.MP_AB = randn(posLength, 1);
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
    cp = max(0, cp,'includenan'); 
end

function [peakIndices, peakValues] = exclusionZonePeaks(profile, mm, exclusionLength, K, distanceThreshold)
%%% Returns the indices of the largest values
%%%   sorted by lowest values
%%%   using excluson zone around matches
%%%   K (optional): return Top-K min indices. Default return all possible indices, 
%%%   WARNING: may not
    numSubsequences = length(profile) - mm + 1;
    if nargin < 3
        exclusionLength = ceil(mm/2);
    end
    if nargin < 4
        K = ceil(2*numSubsequences/exclusionLength);
    end
    if nargin < 5
        distanceThreshold = 0;
    end
    
    
    maxNumPeaks = min(K, ceil(2*numSubsequences/exclusionLength));
    
    
    A = zeros(2,length(profile));
    A(1,:) = 1:length(profile);
    A(2,:) = (profile)';
    A = A';
    B = flip(sortrows(A,2),1);
    
    exclusionZone = zeros(1,length(profile));
    
    peakIndices = nan(maxNumPeaks,1);
    peakValues = nan(maxNumPeaks,1);
    index = 1;
    Kindex = 1; %for debugging
    while index <= maxNumPeaks && Kindex < size(B,1)
        trialIndex = B(Kindex,1);
        if profile(trialIndex) < distanceThreshold
            break;
        end
        if exclusionZone(trialIndex) == 0 && ~isnan(profile(trialIndex))
            peakIndices(index) = trialIndex;
            peakValues(index) = profile(trialIndex);
            index = index + 1;
            leftOffset = exclusionLength;
            rightOffset = exclusionLength;
            exclusionZone(max(1,trialIndex - leftOffset):min(length(profile), trialIndex + rightOffset)) = 1;
        end
        Kindex = Kindex + 1;
    end
    
    peakIndices = peakIndices(~isnan(peakIndices));
    peakValues = peakValues(~isnan(peakIndices));
    
    [peakIndices, sortIndices] = sort(peakIndices);
    peakValues = peakValues(sortIndices);
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