%%%https://www.mathworks.com/matlabcentral/fileexchange/29544-figure-to-play-and-analyze-videos-with-custom-plots-on-top
%%%https://blogs.mathworks.com/pick/2010/12/10/video-player-for-your-frame-based-processing/
% type redraw
function redraw(frame, ssrs, outputVideoHandle)
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
    currIndex = frame + startIndex - 1;
    videoFrames{ii} = ssrs{ii}{2}.read(currIndex);
end


compositeImage = [];
padWidth = 10;
for ii = 1:numSS
    compositeImage = [compositeImage, padarray(videoFrames{ii}, [padWidth,padWidth], nan)];
end

blackHeight = ceil(0.5*size(compositeImage,1));
blackWidth = size(compositeImage,2);
blackBottom = zeros(blackHeight, blackWidth, 3);
compositeImage = cat(1, compositeImage, blackBottom);

% Display
image(compositeImage); 
axis image off;

hold on;
pastWidth = 0;
for ii = 1:numSS
    plotDims = [pastWidth + 1, 1+ssrs{ii}{2}.Height, ssrs{ii}{2}.Width, ceil(0.5*ssrs{ii}{2}.Height)];
    subsequenceColor = ssrs{ii}{4};
    overlayPlot(ssrs{ii}{1}, plotDims, frame, subsequenceColor);
    pastWidth = pastWidth + ssrs{ii}{2}.Width;
end
hold off;

if nargin 
outputFrame = getframe(gcf);
writeVideo(outputVideoHandle, outputFrame);
end

function overlayPlot(sswc, plotDims, frame, subsequenceColor)
    %%%[originX, originY, width, height]
    originX = plotDims(1);
    originY = plotDims(2);
    plotWidth = plotDims(3);
    plotHeight = plotDims(4);

    wcInterpX = originX + linspace(1, plotWidth, sswc.length)-1;
    wcInterpY = interp1([sswc.tsMin,sswc.tsMax],[originY+plotHeight-1, originY],sswc.subsequenceWithContext);
    wocInterpX = interp1(1:sswc.length, wcInterpX, sswc.subsequenceStartIndex:sswc.subsequenceEndIndex)-1;
    wocInterpY = interp1([sswc.tsMin, sswc.tsMax],[originY+plotHeight-1, originY],sswc.subsequenceWithoutContext);

    wcCursorX = interp1([1,sswc.length],[originX,originX+plotWidth-1],frame);
    wcCursorY = interp1([sswc.tsMin, sswc.tsMax], [originY+plotHeight-1, originY], [sswc.subsequenceWithContext(frame)]);

    hold on;
    plot(wcInterpX, wcInterpY, 'Color', [0.7, 0.7, 0.7]);
    plot(wocInterpX, wocInterpY,'Color',subsequenceColor,'LineWidth',3);
    plot(wcCursorX, wcCursorY,'r|','MarkerSize',35,'LineWidth',3);
    
end