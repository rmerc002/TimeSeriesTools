%%%https://www.mathworks.com/matlabcentral/fileexchange/29544-figure-to-play-and-analyze-videos-with-custom-plots-on-top
%%%https://blogs.mathworks.com/pick/2010/12/10/video-player-for-your-frame-based-processing/
% type redraw
function redraw(frame, ssrs, featureAnnotations, outputVideoHandle)
% REDRAW  Process a particular frame of the video
%   REDRAW(FRAME, VIDOBJ)
%       frame  - frame number to process
%       vidObj - VideoReader object

%%% ssrs: Each SubSequence Result Struct is ordered as follows:
%%%     1: sswc
%%%     2: vid
%%%     3: startIndexVideo
%%%     4: subsequenceColor
%%%     5: videoResults

numSS = length(ssrs);
% Read frame
videoFrames = {};
for ii = 1:numSS
    startIndex = ssrs{ii}{3};
    timeIndex = frame + startIndex - 1;
    timeIndex = max(1, timeIndex);
    timeIndex = min(ssrs{ii}{2}.NumFrames, timeIndex);
    videoFrames{ii} = ssrs{ii}{2}.read(timeIndex);
end

%%%Add black area below video to display Contrast Profile plots and text
%%%Also add white border between video
compositeImage = [];
padWidth = 10;
for ii = 1:numSS
    blackHeight = ceil(0.5*size(videoFrames{ii},1) + 0.5*size(videoFrames{ii},1));
    blackWidth = size(videoFrames{ii},2);
    blackBottom = zeros(blackHeight, blackWidth, 3);
    partitionImage = cat(1, videoFrames{ii}, blackBottom);

    %%%white boundary between video
    if ii > 1
        compositeImage = [compositeImage, 255*ones(size(partitionImage,1),padWidth,3)];
    end
    compositeImage = [compositeImage,partitionImage];
end

partitionWidth = size(videoFrames{ii},2) + padWidth;


% Display
image(compositeImage); 
axis image off;

hold on;


plato = ssrs{1}{5}{1}.subsequenceWithoutContext;
platoNN = ssrs{2}{5}{1}.subsequenceWithoutContext;
negNN = ssrs{3}{5}{1}.subsequenceWithoutContext;

distVals = zeros(3,1);
distVals(2) = norm(zscore(plato)-zscore(platoNN))/sqrt(2*length(plato)); %%%Pos Euclidean Dist
distVals(3) = norm(zscore(plato)-zscore(negNN))/sqrt(2*length(plato));%%%Neg Euclidean Dist
distVals(1) = max(0, distVals(3) - distVals(2)); %%%Contrast

distStrings = ["Contrast", "T-Positive Euclidean Distance", "T-Negative Euclidean Distance"];
%%% Text Info
%%% Title
%%% Limb string and feature index
%%% Current behavior Time and index
%%% FileID
titles = ["Most Contrasting Positive Behavior"; "Closest Positive Match"; "Closest Negative Match"];
textInfo = cell(numSS,1);
for ii = 1:numSS
    limbAnn = featureAnnotations{3}(ssrs{ii}{5}{4});
    coordAnn = featureAnnotations{4}(ssrs{ii}{5}{4});
    line1 = sprintf("Limb: %s(%s)    FeatureID: %d",limbAnn, coordAnn, ssrs{ii}{5}{4});
    timeIndex = ssrs{ii}{5}{1}.tsStartIndexContext + frame - 1;
    timeStr = concatenatedIndexToTime(timeIndex, ssrs{ii}{2}.FrameRate);
    line2 = sprintf("Video time: %s    Video frame: %d",timeStr, timeIndex);
    line3 = sprintf("File ID: %s", ssrs{ii}{5}{2});
    line4 = sprintf("%s: %.2f", distStrings(ii), distVals(ii));
    textInfo{ii} = {titles(ii); line1; line2; line3; line4};
end

% textInfo = {'Title'};
pastWidth = 0;

if frame > 1
    for ii = 1:numSS
        plotDims = [pastWidth + 1, 1+ssrs{ii}{2}.Height, ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height), ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height)];
        subsequenceColor = ssrs{ii}{4};
        subsequencePlot(ssrs{ii}{1}, plotDims, frame, subsequenceColor, textInfo{ii});
        pastWidth = pastWidth + partitionWidth;
    end
elseif numSS > 1
    ii = 1;
    plotDims = [pastWidth + 1, 1+ssrs{ii}{2}.Height, ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height), ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height)];
    subsequenceColor = ssrs{ii}{4};
    subsequencePlot(ssrs{ii}{1}, plotDims, frame, subsequenceColor, textInfo{ii});
    pastWidth = pastWidth + partitionWidth;

    for ii = 2:numSS
        plotDims = [pastWidth + 1, 1+ssrs{ii}{2}.Height, ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height), ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height)];
        subsequenceColors = [ssrs{1}{4}; ssrs{ii}{4}];
        overlayPlot(ssrs{1}{1}, ssrs{ii}{1}, plotDims, frame, subsequenceColors, textInfo{ii});
        pastWidth = pastWidth + partitionWidth;
    end
end
hold off;

if nargin == 4
    outputFrame = getframe(gcf);
    writeVideo(outputVideoHandle, outputFrame);
end
end

function subsequencePlot(sswc, plotDims, frame, subsequenceColor, textInfo)
    %%%[originX, originY, width, height]
    originX = plotDims(1);
    originY = plotDims(2);
    plotWidth = plotDims(3);
    plotHeight = plotDims(4);
    textWidth = plotDims(5);
    textHeight = plotDims(6);
    textOriginX = originX;
    textOriginY = originY + plotHeight + ceil(0.5*textHeight);
    textXOffset = 10;

    wcInterpX = originX + linspace(1, plotWidth, sswc.length)-1;
    wcInterpY = interp1([sswc.tsMin,sswc.tsMax],[originY+plotHeight-1, originY],sswc.subsequenceWithContext);
    wocInterpX = interp1(1:sswc.length, wcInterpX, sswc.subsequenceStartIndex:sswc.subsequenceEndIndex)-1;
    wocInterpY = interp1([sswc.tsMin, sswc.tsMax],[originY+plotHeight-1, originY],sswc.subsequenceWithoutContext);

    wcCursorX = interp1([1,sswc.length],[originX,originX+plotWidth-1],frame);
    wcCursorY = interp1([sswc.tsMin, sswc.tsMax], [originY+plotHeight-1, originY], [sswc.subsequenceWithContext(frame)]);

    hold on;
    plot(wcInterpX, wcInterpY, 'Color', [0.7, 0.7, 0.7]);
    text(textOriginX + textXOffset, textOriginY,textInfo,'Color',[1,1,1],'VerticalAlignment', 'middle');

    plot(wocInterpX, wocInterpY,'Color',subsequenceColor,'LineWidth',3);

    plot(wcCursorX, wcCursorY,'r|','MarkerSize',35,'LineWidth',1);
    
end

function overlayPlot(sswc1, sswc2, plotDims, frame, subsequenceColors, textInfo)
    %%%[originX, originY, width, height]
    originX = plotDims(1);
    originY = plotDims(2);
    plotWidth = plotDims(3);
    plotHeight = plotDims(4);
    textWidth = plotDims(5);
    textHeight = plotDims(6);
    textOriginX = originX;
    textOriginY = originY + plotHeight + ceil(0.5*textHeight);
    textXOffset = 10;

    ss1 = zscore(sswc1.subsequenceWithoutContext');
    ss1 = reshape(ss1, [1,length(ss1)]);
    ss2 = zscore(sswc2.subsequenceWithoutContext');
    ss2 = reshape(ss2,[1,length(ss2)]);

    ssMin = min([ss1,ss2]);
    ssMax = max([ss1, ss2]);

    wocInterpX1 = originX + linspace(1, plotWidth, length(ss1))-1;
    wocInterpY1 = interp1([ssMin, ssMax],[originY+plotHeight-1, originY],ss1);

    wocInterpX2 = wocInterpX1;
    wocInterpY2 = interp1([ssMin, ssMax],[originY+plotHeight-1, originY],ss2);

    wcCursorX = interp1([1,sswc1.length],[originX,originX+plotWidth-1],frame);
    wcCursorY = interp1([ssMin, ssMax], [originY+plotHeight-1, originY], [sswc1.subsequenceWithContext(frame)]);

    hold on;
    plot(wocInterpX1, wocInterpY1,'Color',subsequenceColors(1,:),'LineWidth',3);
    plot(wocInterpX2, wocInterpY2,'Color',subsequenceColors(2,:),'LineWidth',3);

    text(textOriginX + textXOffset, textOriginY,textInfo,'Color',[1,1,1],'VerticalAlignment', 'middle');

%     plot(wcCursorX, wcCursorY,'r|','MarkerSize',35,'LineWidth',1);
    
end