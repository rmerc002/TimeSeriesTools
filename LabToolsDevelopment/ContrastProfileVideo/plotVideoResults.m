function plotVideoResults(videoResults, outputPath)
%%%Expect results package in the form of 
%%% level 1: {plato, platoNN, negNN}
%%% level 2: {subsequenceObj; fileNumber; videoPath; posNegFeatureIndices, tsMinMax}
%%% subsequnceObj properties:
%%%     subsequenceWithContext;
%%%     contextLength;
%%%     subsequenceLength;
%%%     length;
%%%     tsMin;
%%%     tsMax;
%%%     tsStartIndexContext;
%%%     tsStartIndexSubsequence;


%%%There will be a problem if the subsequence found is less than context
%%%length away from either end of the original ts.

%%%Dependencies:
%%% sbusequence.m
%%% redraw.m
%%% videofig

sswcPlato = videoResults{1}{1};
vidPlato = VideoReader(videoResults{1}{3});
startIndexVideoPlato = sswcPlato.tsStartIndexContext;
numFrames = sswcPlato.length;%vid1.NumberOfFrames
subsequenceColor = [129/255, 51/255, 144/255]; %%% purple
platoStruct = {sswcPlato, vidPlato, startIndexVideoPlato, subsequenceColor, videoResults{1}};

sswcPlatoNN = videoResults{2}{1};
vidPlatoNN = VideoReader(videoResults{2}{3});
startIndexVideoPlatoNN = sswcPlatoNN.tsStartIndexContext;
subsequenceColor = [115/255, 170/255, 43/255]; %%% green
platoNNStruct = {sswcPlatoNN, vidPlatoNN, startIndexVideoPlatoNN, subsequenceColor, videoResults{2}};

sswcNegNN = videoResults{3}{1};
vidNegNN = VideoReader(videoResults{3}{3});
startIndexVideoNegNN = sswcNegNN.tsStartIndexContext;
subsequenceColor = [0.73,0.05,0]; %%% red
negNNStruct = {sswcNegNN, vidNegNN, startIndexVideoNegNN, subsequenceColor, videoResults{3}};

subsequenceResultStructs = {platoStruct, platoNNStruct, negNNStruct};

if nargin == 2
    outputVideoHandle = VideoWriter(outputPath, 'MPEG-4');
    outputVideoHandle.FrameRate = 30;
    open(outputVideoHandle);
    videofig(numFrames, @(frm) redraw(frm, subsequenceResultStructs, outputVideoHandle));
else
    videofig(numFrames, @(frm) redraw(frm, subsequenceResultStructs));
    redraw(1, subsequenceResultStructs);
end

fig = gcf;
scale = 200;
% fig.Position = [100,100,ceil(scale*16),ceil(scale*2*9)];
fig.Units = 'pixels';
fig.Position = [760,600,1275,500];

if nargin == 2
    for ii = 1:numFrames
        redraw(ii, subsequenceResultStructs, outputVideoHandle);
    end
    close(outputVideoHandle)
    close(gcf);
end

end